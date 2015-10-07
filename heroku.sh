#!/bin/bash
gunicorn flory:app
python worker.py
