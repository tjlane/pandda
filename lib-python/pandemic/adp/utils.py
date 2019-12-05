
def show_file_dict(self, file_dict, indent=0):
    log = self.log
    s = '  '
    for k, v in file_dict.iteritems():
        if isinstance(v, dict):
            log(s*indent + '> {}'.format(k))
            show_file_dict(self, v, indent+1)
        elif isinstance(v, str):
            log(s*indent + '> {}: {}'.format(k, v))
        else:
            log(s*indent + '> {}'.format(k))
            try:
                for vv in v:
                    log(s*(indent+1)+str(vv))
            except:
                log(s*(indent+1)+str(v))