from google.appengine.ext import webapp

class Guesser(webapp.RequestHandler):
    formstring = '''<form method="post" action="/">
<p>Enter Guess: <input type="tesxt" name="guess"/></p>
<p><input type="submit"></p>
</form>'''

    def get(self):
        self.response.out.write('<p>Good luck!</p>\n')
        self.response.out.write(self.formstring)

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

        self.response.out.write('<p>Guess:'+stguess+'</p>\n')
        self.response.out.write('<p>'+msg+'</p>')
        self.response.out.write(self.formstring)
