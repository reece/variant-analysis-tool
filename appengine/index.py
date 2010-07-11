import wsgiref.handlers
from google.appengine.ext import webapp

from vat import VAT, TabHandler, Redirect

def main():
    application = webapp.WSGIApplication([
		('/tab/.*', TabHandler),
		('/vat.*', VAT),
		('/', Redirect)],
		debug=True)						# print debug info to browser
    wsgiref.handlers.CGIHandler().run(application)

if __name__ == '__main__':
    main()
