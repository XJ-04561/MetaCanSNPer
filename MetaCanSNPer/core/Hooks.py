
import itertools
from typing import Callable, Any, overload
from functools import cache
from threading import Thread, current_thread
from queue import Queue, Empty as EmptyQueueException

from MetaCanSNPer.Globals import Logged, forceHash

class HookBaitingException(Exception): pass

class Hook:

	target : Callable
	args : tuple
	kwargs : dict

	def __init__(self, target, args=(), kwargs={}):
		self.target = target
		self.args = args
		self.kwargs = kwargs
	
	def __call__(self, eventInfo : dict):
		self.target(eventInfo, *self.args, **self.kwargs)
	
	def __hash__(self):
		return forceHash((self.target, self.args, self.kwargs))
	
	def __eq__(self, other):
		return hash(self) == hash(other)

class DummyHooks:

	_eventQueue : Queue[tuple[str,dict]]
	_hooks : dict[str, set[Hook]]
	_worker : Thread
	RUNNING : bool

	def addHook(self, eventType : str, target : Callable, args : tuple=(), kwargs : dict={}) -> Hook:
		return Hook(lambda *args, **kwargs : None, args=args if hasattr(args, "__hash__") else tuple(args), kwargs=kwargs)
	
	def removeHook(self, eventType : str, hook : Hook) -> bool:
		return True

	def trigger(self, eventType : str, eventInfo : dict):
		pass

class Hooks(Logged):
	"""
	eventInfo should look like:
		*name* - This can be any string. But is used together with the identity of the 'instance' to make the trigger precise in its targeting.
		
		*instance* - The instance which triggered the event, or the instance relevant to the event.
		
		*owner* - The class type of the 'instance'.
		
		*value* - The value which the event is reporting.
		
		_ - All other keywords are available for use, but should not be on relied on for identification or value-passing when the 'value' key is available.
	"""

	_eventQueue : Queue[tuple[str,dict]]
	_hooks : dict[str, set[Hook]]
	_worker : Thread
	RUNNING : bool

	def __init__(self):
		self._hooks = {}
		self._eventQueue = Queue()
		self.RUNNING = True
		self._worker = Thread(target=self.mainLoop, daemon=True)
		self._worker.start()
	
	def __del__(self):
		self.RUNNING = False
	@overload
	def addHook(self, eventType : str, hook : Hook) -> Hook: ...
	@overload
	def addHook(self, eventType : str, target : Callable, args : tuple=(), kwargs : dict={}) -> Hook: ...
	def addHook(self, eventType : str, target : Callable, args : tuple=(), kwargs : dict={}) -> Hook:
		
		if isinstance(target, Hook):
			hook = target
		else:
			args = args if hasattr(args, "__hash__") else tuple(args)
			self.LOG.debug(f"Adding Hook to {self}. Hook has: {eventType=}, {target=}, {args=}, {kwargs=}")
			
			hook = Hook(target, args, kwargs)
		
		if eventType not in self._hooks:
			self._hooks[eventType] = set()
		self._hooks[eventType].add(hook)
		
		return hook
	
	def removeHook(self, eventType : str, hook : Hook) -> bool:
		"""Removes all occurances of the hook in the list of hooks associated with the eventType"""
		# Essentially a while True: but limited to at least iterations as long as the hooks list.
		if hook in self._hooks.get(eventType, set()):
			self._hooks[eventType].remove(hook)
			return True
		else:
			return False
	
	def trigger(self, eventType : str|tuple, eventInfo : dict):
		
		if isinstance(eventType, str):
			self.LOG.debug(f"Event triggered: {eventType=}, {eventInfo=}")
			self._eventQueue.put((eventType, eventInfo))
		elif isinstance(eventType, tuple):
			for name in eventType:
				self.LOG.debug(f"Event triggered: {name=}, {eventInfo=}")
				self._eventQueue.put((name, eventInfo))
		else:
			self.LOG.exception(TypeError(f"eventType must be either a str or a tuple of str. not {eventType!r}\n{eventInfo=}"))
			raise TypeError(f"eventType must be either a str or a tuple of str. not {eventType!r}\n{eventInfo=}")
	
	def mainLoop(self):
		
		while self.RUNNING:
			try:
				eventType, eventInfo = self._eventQueue.get(timeout=2)
				if self._hooks.get(eventType): 
					for hook in itertools.takewhile(lambda x:self.RUNNING, self._hooks.get(eventType, [])):
						try:
							hook(eventInfo)
						except Exception as e:
							e.add_note(f"This occurred while calling the hook {hook!r} tied to {eventType=} with {eventInfo=}")
							self.LOG.exception(e)
				else:
					self.LOG.warning(f"{eventType=} triggered, but no hooks registered in {self!r}")
			except EmptyQueueException:
				pass
			except Exception as e:
				e.add_note(f"This exception occurred in hooks thread '{getattr(current_thread(), 'name', 'N/A')}'")
				self.LOG.exception(e)
		self.LOG.info(f"Hooks thread {getattr(current_thread(), 'name', 'N/A')} stopped running")


class Bait(Logged):
	"""

	"""

	name : str = None
	eventType : str|tuple[str] = None
	eventInfo : dict = None
	owner : type = None
	attributeName : str = None
	_property : property = None

	def __init__(self, *, name : str=None, eventType : str=None, eventInfo : dict=None):
		
		self.eventType = eventType or self.eventType
		self.eventInfo = eventInfo or self.eventInfo or {}
		self.name = name or self.name
		self.eventInfo.setdefault("name", self.name)
	
	def __class_getitem__(cls, name : str):
		
		return cls(name=name)
	
	def __getitem__(self, eventTypes : str|tuple[str]):
		
		return type(self)(name=self.name, eventType=eventTypes)
	
	def __call__(self, fget=None, fset=None, fdel=None, doc=None):
		
		if isinstance(fget, property):
			self._property = fget
		else:
			self._property = property(fget=fget, fset=fset, fdel=fdel, doc=doc)
	
	def __get__(self, instance):
		
		if self.name in instance.__dict__:
			return instance.__dict__[self.name]
		else:
			raise AttributeError(f"{type(instance).__name__!r} object has no attribute {self.name!r}")

	def __set__(self, instance, value):
		
		instance.__dict__[self.name] = value
		getattr(instance, "hooks").trigger(self.eventType, self.eventInfo | {"value" : value, "instance" : instance})
	
	def __set_name__(self, owner : type, name : str):
		
		if not hasattr(owner, "hooks") and "hooks" not in owner.__annotations__ and "hooks" not in owner.__init__.__code__.co_varnames[:owner.__init__.__code__.co_argcount]:
			raise HookBaitingException("Can't lay bait on an attribute belonging to a class that does not have a 'hooks' attribute. ('hooks' is determined through class attribute lookup, class annotation lookup, and __init__ function arguments names)")
		self.eventType = self.eventType or f"{owner.__name__}{name.capitalize()}"
		self.eventInfo["owner"] = self.owner = owner
		self.eventInfo["attribute"] = self.attributeName = name
		self.eventInfo.setdefault("name", name)

	
	def __delete__(self, instance):
		
		del instance.__dict__[self.name]

	def __repr__(self):
		
		return f"<Baited {type(self._property).__name__ if self._property else 'Attribute'} '{self.owner.__name__}.{self.attributeName}' on {self.eventType!r} with {self.eventInfo!r}>"
	
	def setter(self, func):
		
		if self._property:
			self._property.setter(func)
		else:
			self._property = property(fset=func)
		
	def deleter(self, func):
		
		if self._property:
			self._property.deleter(func)
		else:
			self._property = property(fdel=func)

GlobalHooks = Hooks()