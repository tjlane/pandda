
def select_factors_and_phases_for_map(cols, maptype='2FOFC'):
    """Takes the Column Headings and returns a dictionary for creating a map type"""

    # Check for Valid Map
    if maptype not in ['2FOFC', 'FOFC']:
        raise Exception('MapType not in allowed map types: {!s} not in {!s}'.format(maptype,MAP_TYPES))
    # Dictionary to contain headings
    lab_dict = {'F':'','PHI':''}
    # Match
    if maptype == '2FOFC':
        if ('2FOFCWT' in cols) and ('PH2FOFCWT' in cols):
            lab_dict['F'] = '2FOFCWT'
            lab_dict['P'] = 'PH2FOFCWT'
        elif ('FWT' in cols) and ('PHWT' in cols):
            lab_dict['F'] = 'FWT'
            lab_dict['P'] = 'PHWT'
        elif ('FWT' in cols) and ('PHFWT' in cols):
            lab_dict['F'] = 'FWT'
            lab_dict['P'] = 'PHFWT'
        else:
            raise ValueError('Failed to select STRUCTURE FACTOR Columns for Map Type {!s}. \nColumns : {!s}'.format(maptype,cols))
    elif maptype == 'FOFC':
        if ('FOFCWT' in cols) and ('PHFOFCWT' in cols):
            lab_dict['F'] = 'FOFCWT'
            lab_dict['P'] = 'PHFOFCWT'
        elif ('DELFWT' in cols) and ('DELPHWT' in cols):
            lab_dict['F'] = 'DELFWT'
            lab_dict['P'] = 'DELPHWT'
        elif ('DELFWT' in cols) and ('PHDELWT' in cols):
            lab_dict['F'] = 'DELFWT'
            lab_dict['P'] = 'PHDELWT'
        else:
            raise ValueError('Failed to select PHASE Columns: {!s}'.format(cols))
    else:
        raise NotImplementedError('MapType does not have functionality programmed: {!s}'.format(maptype))
    # Return
    return lab_dict

