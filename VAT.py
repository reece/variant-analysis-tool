import os
import logging
import wsgiref.handlers
from google.appengine.ext import webapp
from google.appengine.ext.webapp import template
from Bio import Entrez

class VAT(webapp.RequestHandler):
    def get(self):
        path = self.request.path
        temp = os.path.join(
            os.path.dirname(__file__),
            'templates/VAT.html')
        self.response.out.write( template.render(
                temp,
                { 
                    'path': path,
                    'pagetitle': 'Variant Analysis Tool' 
                    }
                ))

    def post(self):
        stguess = self.request.get('guess')
        try:
            guess = int(stguess)
        except:
            guess = -1

        answer = 42
        if guess == answer:
            msg = 'Congratulations!'
        elif guess < 0:
            msg = 'Please provide a number guess'
        elif guess < answer:
            msg = 'Your guess is too low'
        else:
            msg = 'Your guess is too high'

        self.response.out.write('<h1>VAT</h1>>\n')
        self.response.out.write('<p>Guess:'+stguess+'</p>\n')
        self.response.out.write('<p>'+msg+'</p>')
        self.response.out.write(self.formstring)
