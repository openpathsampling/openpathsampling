from celery import Celery
from celery.signals import worker_process_init
import logging

logger = logging.getLogger(__name__)


app = Celery('tasks', backend='redis://localhost:6379/0', broker='amqp://')
app.config_from_object('celeryconfig')


# This will force each child process to have a sep

def _redo_uuid(**kwargs):
    import openpathsampling.netcdfplus as npl
    # logger.info('OLD UUID `%s`' % npl.StorableObject.INSTANCE_UUID)
    npl.StorableObject.initialize_uuid()
    logger.info('NEW UUID `%s`' % npl.StorableObject.INSTANCE_UUID)

# need to use registered instance for sender argument.
worker_process_init.connect(_redo_uuid)


@app.task
def add(x, y):
    return x + y


@app.task
def generate(engine, template, ensemble):
    engine.initialize('CPU')
    traj = engine.generate(template, ensemble.can_append)
    return traj


@app.task
def uuid(obj):
    return obj.__uuid__


@app.task
def identity(obj):
    return obj


@app.task
def attr(obj, attr):
    return getattr(obj, attr)


@app.task
def ls(obj):
    return dir(obj)


@app.task
def temp(obj):
    return obj._lazy.values()[0].__dict__


@app.task
def exev(template, s1, s2):
    exec(s1)
    return eval(s2)

