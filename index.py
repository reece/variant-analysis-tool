import wsgiref.handlers
from google.appengine.ext import webapp

from VAT import VAT


def main():
    application = webapp.WSGIApplication([
            ('/.*', VAT)],
            debug=True)
    wsgiref.handlers.CGIHandler().run(application)


if __name__ == '__main__':
    main()
