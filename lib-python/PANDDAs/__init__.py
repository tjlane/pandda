

def welcome(current_user):

    try:
        from greetings import get_greeting
        print get_greeting(username=current_user)
        return
    except:
        pass

    print 'Hi {!s}. Welcome to PANDDAs.'.format(current_user.upper())
