from celery import Celery
from celery.signals import worker_process_init, celeryd_init
import logging

from sshtunnel import SSHTunnelForwarder

from celery.signals import worker_init, worker_shutdown

logger = logging.getLogger(__name__)

db_server = 'shark.imp.fu-berlin.de'
redis_db_port = 6379

redis_server = (db_server, 22)
node_remote = ('localhost', redis_db_port)

node_port = 6383  # can be any port. Needs to be the same as in redis worker

redis_server_user = 'jprinz'
keyfile = '/Users/jan-hendrikprinz/.ssh/known_hosts'

ssh_password = open('pw').read()

server_str = 'redis://localhost:%d/0' % node_port

app = Celery('tasks', backend=server_str, broker=server_str)
app.config_from_object('celeryconfig')

tunnel_server = SSHTunnelForwarder(
    db_server,
    # ssh_host_key=known_hosts_line,
    ssh_username=redis_server_user,
    ssh_password=ssh_password,
    local_bind_address=('127.0.0.1', node_port),
    remote_bind_address=('127.0.0.1', redis_db_port)
)


# This will force each child process to have a seperate INSTANCE_UUID
# only necessary if you fork your processes in which case python
# will reuse exising variables that _seem_ untouched and the ops package
# import will not be executed again
# on remote machine this should actually not be necessary

def _redo_uuid(**kwargs):
    import openpathsampling.netcdfplus as npl
    # logger.info('OLD UUID `%s`' % npl.StorableObject.INSTANCE_UUID)
    npl.StorableObject.initialize_uuid()
    logger.info('NEW UUID `%s`' % npl.StorableObject.get_instance_uuid())

# need to use registered instance for sender argument.
worker_process_init.connect(_redo_uuid)


def _open_tunnel(**kwargs):
    logger.info('RUNNING TUNNEL')
    tunnel_server.start()
    logger.info('tunnel established at port %s' % str(tunnel_server.local_bind_address))


def _close_tunnel(**kwargs):
    logger.info('Shutting down tunnel')
    tunnel_server.stop()
    logger.info('Tunnel removed')


celeryd_init.connect(_open_tunnel)
worker_shutdown.connect(_close_tunnel)


@app.task(name='openpathsampling.engine.celery.tasks.generate')
def generate(engine, template, ensemble, init_args=None, init_kwargs=None):
    if init_args is None:
        init_args = []
    if init_kwargs is None:
        init_kwargs = {}

    engine.initialize(*init_args, **init_kwargs)
    traj = engine.generate(template, ensemble.can_append)
    return traj


@app.task(name='openpathsampling.engine.celery.tasks.list_cache')
def list_cache():
    """
    Return a list of UUIDs containing all objects present in the cache

    This is useful for checking which objects should actually be
    transferred. This is not save yet, becuase in the meantime the
    cache could loose some objects in which case you will not transmit these.
    This needs a function to lock the cache in that time and release after
    transmission


    Returns
    -------
    list of UUID
        the list of UUIDs of the objects that do not beed to be transmitted
    """
    return []

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

