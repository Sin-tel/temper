from contextlib import contextmanager

USE_TIMEOUT = True


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    if USE_TIMEOUT:
        import signal

        def signal_handler(signum, frame):
            raise TimeoutException("Timed out!")

        if not hasattr(signal, "SIGALRM"):
            print("SIGALRM not available!")
            return

        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)
    else:
        pass
