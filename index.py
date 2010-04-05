import wsgiref.handlers
from google.appengine.ext import webapp
from Guesser import Guesser


def main():
    application = webapp.WSGIApplication([
            ('/.*', Guesser)],
            debug=True)
    wsgiref.handlers.CGIHandler().run(application)


if __name__ == '__main__':
    main()
