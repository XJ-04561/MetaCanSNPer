from MetaCanSNPer.core.Hooks import Hooks, Hook



def test_hooks():

	h = Hook(target=lambda *args, **kwargs: print("Hej!"))

	hooks = Hooks()

	hooks.addHook("testHook", h)

	hooks.trigger("testHook", {})

	assert True == hooks.removeHook("testHook", h)