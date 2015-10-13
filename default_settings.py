"""Default Flask app configuration.

Intended to be loaded with::

    app.config.from_object(rq_dashboard.default_settings)

and then overridden as required.

"""
#: If set the REDIS_URL takes precedence over REDIS_HOST, REDIS_PORT, etc
REDIS_URL = redis://redistogo:d302ccb0f75a68e30cead681470e6588@cowfish.redistogo.com:9273/

REDIS_HOST = 'localhost'
REDIS_PORT = 6379
REDIS_PASSWORD = None
REDIS_DB = 0

RQ_POLL_INTERVAL = 2500  #: Web interface poll period for updates in ms
DEBUG = True
