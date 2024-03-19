
from typing import Callable
from threading import Thread, Condition

from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger(__name__)

class Hook:

    target : Callable
    args : list
    kwargs : dict

    def __init__(self, target, args=[], kwargs={}):
        self.target = target
        self.args = args
        self.kwargs = kwargs
    
    def __call__(self, eventInfo : dict):
        self.target(eventInfo, *self.args, **self.kwargs)

class Hooks:

    _eventQueue : list[tuple[str,dict]]
    _hooks : dict[str, list[Hook]]
    _worker : Thread
    RUNNING : bool

    def __init__(self):
        self._hooks = {}
        self._eventQueue = []
        self.RUNNING = True
        self._worker = Thread(target=self.mainLoop, daemon=True)
        self._worker.start()
    
    def __del__(self):
        self.RUNNING = False
    
    def addHook(self, eventType : str, target : Callable, args : list=[], kwargs : dict={}) -> Hook:
        
        LOGGER.debug(f"Adding Hook to {self}. Hook has: {eventType=}, {target=}, {args=}, {kwargs=}")
        if eventType not in self._hooks:
            self._hooks[eventType] = []
        
        hook = Hook(target, args, kwargs)
        self._hooks[eventType].append(hook)

        return hook
    
    def removeHook(self, eventType : str, hook : Hook):
        """Removes all occurances of the hook in the list of hooks associated with the eventType"""
        # Essentially a while True: but limited to at least iterations as long as the hooks list.
        ret = False
        for _ in range(len(self._hooks[eventType])):
            try:
                self._hooks[eventType].remove(hook)
                ret = True
            except ValueError:
                return ret

    def trigger(self, eventType : str, eventInfo : dict):
        LOGGER.debug(f"Event triggered: {eventType=}, {eventInfo=}")
        self._eventQueue.append((eventType, eventInfo))
    
    def mainLoop(self):
        cond = Condition()
        while self.RUNNING:
            cond.acquire(False)
            cond.wait_for(lambda :len(self._eventQueue)>0, timeout=1)
            if len(self._eventQueue)>0:
                eventType, eventInfo = self._eventQueue.pop(0)
                if eventType not in self._hooks or self._hooks.get(eventType) == []: LOGGER.warning(f"{eventType=} triggered, but no hooks registered in {self!r}")
                for hook in self._hooks.get(eventType, []):
                    if self.RUNNING: hook(eventInfo)
