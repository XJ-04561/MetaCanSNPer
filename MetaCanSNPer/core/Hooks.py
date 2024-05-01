
from typing import Callable, Any, overload
from threading import Thread, Condition
from types import GenericAlias

from MetaCanSNPer.core.LogKeeper import createLogger

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

class DummyHooks:

    _eventQueue : list[tuple[str,dict]]
    _hooks : dict[str, list[Hook]]
    _worker : Thread
    RUNNING : bool

    def addHook(self, *args, **kwargs) -> Hook:
        return Hook(lambda *args, **kwargs : None)
    
    def removeHook(self, eventType : str, hook : Hook) -> bool:
        return True

    def trigger(self, eventType : str, eventInfo : dict):
        pass

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
    
    def removeHook(self, eventType : str, hook : Hook) -> bool:
        """Removes all occurances of the hook in the list of hooks associated with the eventType"""
        # Essentially a while True: but limited to at least iterations as long as the hooks list.
        ret = False
        if eventType not in self._hooks:
            return ret
        for _ in range(len(self._hooks[eventType])):
            try:
                self._hooks[eventType].remove(hook)
                ret = True
            except ValueError:
                return ret
        return ret

    def trigger(self, eventType : str|tuple, eventInfo : dict):
        if isinstance(eventType, str):
            LOGGER.debug(f"Event triggered: {eventType=}, {eventInfo=}")
            self._eventQueue.append((eventType, eventInfo))
        else:
            for name in eventType:
                LOGGER.debug(f"Event triggered: {name=}, {eventInfo=}")
                self._eventQueue.append((name, eventInfo))
    
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

class Trigger:

    owner : type
    instance : int
    name : str
    eventType : str
    eventInfo : dict
    nameAttr : str = ""
    values : dict[int,Any]
    def __init__(self, eventType : str=None, nameAttr : str=None, eventInfo : dict={}):
        self.eventType = eventType
        self.eventInfo = eventInfo
        self.nameAttr = nameAttr or self.nameAttr
        self.values = {}
    
    def __get__(self, instance, owner=None):
        return self.values[id(instance)]

    def __set__(self, instance, value):
        self.values[id(instance)] = value
        instance.hooks.trigger(self.eventType, self.eventInfo | {self.name : value, "instance" : instance, "name" : getattr(instance, self.self.nameAttr, None)})
    
    def __set_name__(self, owner : type, name : str):
        self.eventType = self.eventType or f"{owner.__name__}{name.capitalize()}"
        self.owner = owner.__name__
        self.name = name
        if owner.__annotations__.get("hooks") is not Hooks:
            raise TypeError(f"{owner!r} has no properly annotated attribute 'hooks'")

    def __repr__(self):
        return f"<Trigger-Property {self.owner}.{self.name!r} on {self.eventType!r}>"

def urlretrieveReportHook(category, hooks : Hooks, name : str, steps : int=-1):

    if steps > 0:
        # info : (blocks, blockSize, totalSize)
        while (info := (yield)) and info[0] * info[1] < info[2]:
            hooks.trigger(f"{category}Progress", {"name" : name, "progress" : min(1.0, info[0] * info[1] / info[2])})
        else:
            hooks.trigger(f"{category}Progress", {"name" : name, "progress" : 1.0})
    else:
        # info : (blocks, blockSize, totalSize)
        pos = 1
        while (info := (yield)) and info[0] * info[1] < info[2]:
            if pos <= info[0] * info[1] / info[2]:
                hooks.trigger(f"{category}", {"name" : name, "progress" : min(1.0, info[0] * info[1] / info[2])})
                pos += 1
        else:
            if pos <= info[0] * info[1] / info[2]:
                hooks.trigger(f"{category}", {"name" : name, "progress" : 1.0})
                pos += 1