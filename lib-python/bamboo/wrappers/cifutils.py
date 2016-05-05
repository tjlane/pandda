import os, shutil, tempfile

from bamboo.common.command import CommandManager

def merge_cif_libraries(incifs, outcif):
    """Take a list of cifs and merge into one cif"""

    assert isinstance(incifs, list), "'incifs' is not a list!"
    assert len(incifs) > 1, "'incifs' must be two or more files!"

    current = incifs.pop(0)

    to_be_deleted = []

    for additional in incifs:

        # Create a file handle and path for the output
        temp_handle, temp_path = tempfile.mkstemp(suffix='.lib', prefix='libcheck_')
        to_be_deleted.append(temp_path)

        # Merge current and additional to temp_path
        LIBCHK = CommandManager('libcheck')
        LIBCHK.add_standard_input('_DOC N','_FILE_L {!s}'.format(current),'_FILE_L2 {!s}'.format(additional),'_FILE_O {!s}'.format(temp_path.replace('.lib','')),'_END')
        LIBCHK.run()

        current = temp_path

    shutil.copy(current, outcif)

    for file in to_be_deleted:
        os.remove(file)

    assert os.path.exists(outcif), 'OUTPUT CIF DOES NOT EXIST! {!s}'.format(outcif)

    return outcif

