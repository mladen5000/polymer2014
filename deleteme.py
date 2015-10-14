from rq import Queue


queue_name = 'lol'
q2 = Queue(queue_name)
q2.get_jobs(1,5)

