
def welcome(current_user):

    try:
        from greetings import get_greeting
        print get_greeting(username=current_user)
        return
    except SystemExit:
        raise
    except:
        pass

    print 'Hi {!s}. Welcome to Pandda.'.format(current_user.upper())