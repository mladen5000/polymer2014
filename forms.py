#!/usr/bin/env python
from flask.ext.wtf import Form

from wtforms import StringField, BooleanField, TextAreaField,SubmitField, validators
from wtforms.validators import DataRequired

class ContactForm(Form):
  name = StringField("Name", [validators.DataRequired()] )
  email = StringField("Email", [validators.DataRequired(), validators.Email()] )
  subject = StringField("Subject", [validators.DataRequired()] )
  message = TextAreaField("Message", [validators.DataRequired()] )
  submit = SubmitField("Send", [validators.DataRequired()] )
