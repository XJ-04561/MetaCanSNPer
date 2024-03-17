
from typing import Callable
from threading import Thread, Condition

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
    _hooks : dict
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
    
    def addHook(self, eventType, target, args=[], kwargs={}):
        if eventType not in self._hooks:
            self._hooks[eventType] = []
        
        self._hooks[eventType].append(Hook(target, args, kwargs))
    
    def trigger(self, eventType, eventInfo):
        self._eventQueue.append((eventType, eventInfo))
    
    def mainLoop(self):
        cond = Condition()
        while self.RUNNING:
            cond.acquire(False)
            cond.wait_for(lambda :len(self._eventQueue)>0, timeout=1)
            if len(self._eventQueue)>0:
                eventType, eventInfo = self._eventQueue.pop(0)
                for hook in self._hooks.get(eventType, []):
                    if self.RUNNING: hook(eventInfo)
