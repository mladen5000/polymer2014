#!/bin/bash
gunicorn flory:app --daemon
python worker.py
