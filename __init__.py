from __future__ import unicode_literals

from math import *

import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins

from flask import Flask, request, make_response, render_template,jsonify
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import json
import subprocess

#Used for Self Consistent Field Theory
import subprocess
from scft import *

#Used to Debug
from flask_debugtoolbar import DebugToolbarExtension

#Task queueing, redis
from rq import Queue
from rq.job import Job
from worker import conn
from rq_dashboard import RQDashboard

#For talking to chiSQL
import json
import urllib
