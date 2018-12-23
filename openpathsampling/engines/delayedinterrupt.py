import signal
import logging


# class based on: http://stackoverflow.com/a/21919644/487556
class DelayedInterrupt(object):
    def __init__(self, signals=None):
        if signals is None:
            signals = [signal.SIGINT]
        if not isinstance(signals, list) and not isinstance(signals, tuple):
            signals = [signals]
        self.sigs = signals

    def __enter__(self):
        self.signal_received = {}
        self.old_handlers = {}
        for sig in self.sigs:
            self.signal_received[sig] = False
            self.old_handlers[sig] = signal.getsignal(sig)

            def handler(s, frame):
                self.signal_received[sig] = (s, frame)
                # Note: in Python 3.5, you can use signal.Signals(sig).name
                logging.info(
                    'Signal %s received. Delaying KeyboardInterrupt.' % sig)

            self.old_handlers[sig] = signal.getsignal(sig)
            signal.signal(sig, handler)

    def __exit__(self, type, value, traceback):
        for sig in self.sigs:
            signal.signal(sig, self.old_handlers[sig])
            if self.signal_received[sig] and self.old_handlers[sig]:
                self.old_handlers[sig](*self.signal_received[sig])

class EmptyContext(object):
    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        pass
