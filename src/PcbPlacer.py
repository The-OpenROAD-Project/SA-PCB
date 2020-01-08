# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.10
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_PcbPlacer', [dirname(__file__)])
        except ImportError:
            import _PcbPlacer
            return _PcbPlacer
        if fp is not None:
            try:
                _mod = imp.load_module('_PcbPlacer', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _PcbPlacer = swig_import_helper()
    del swig_import_helper
else:
    import _PcbPlacer
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class GridBasedPlacer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, GridBasedPlacer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, GridBasedPlacer, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _PcbPlacer.new_GridBasedPlacer(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _PcbPlacer.delete_GridBasedPlacer
    __del__ = lambda self : None;
    def test_placer_flow(self) -> "kicadPcbDataBase &" : return _PcbPlacer.GridBasedPlacer_test_placer_flow(self)
    def getDb(self) -> "kicadPcbDataBase &" : return _PcbPlacer.GridBasedPlacer_getDb(self)
    __swig_setmethods__["mDb"] = _PcbPlacer.GridBasedPlacer_mDb_set
    __swig_getmethods__["mDb"] = _PcbPlacer.GridBasedPlacer_mDb_get
    if _newclass:mDb = _swig_property(_PcbPlacer.GridBasedPlacer_mDb_get, _PcbPlacer.GridBasedPlacer_mDb_set)
    def get_total_cost(self) -> "double" : return _PcbPlacer.GridBasedPlacer_get_total_cost(self)
    def get_wirelength_cost(self) -> "double" : return _PcbPlacer.GridBasedPlacer_get_wirelength_cost(self)
    def get_overlap_cost(self) -> "double" : return _PcbPlacer.GridBasedPlacer_get_overlap_cost(self)
    def get_temperature(self) -> "double" : return _PcbPlacer.GridBasedPlacer_get_temperature(self)
    def set_overlap_weight(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_overlap_weight(self, *args)
    def set_wirelength_weight(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_wirelength_weight(self, *args)
    def set_two_sided(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_two_sided(self, *args)
    def set_initial_move_radius(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_initial_move_radius(self, *args)
    def set_rtree(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_rtree(self, *args)
    def set_lam(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_lam(self, *args)
    def set_lamtemp_update(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_lamtemp_update(self, *args)
    def set_num_iterations(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_num_iterations(self, *args)
    def set_iterations_moves(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_iterations_moves(self, *args)
    def set_initial_temperature(self, *args) -> "void" : return _PcbPlacer.GridBasedPlacer_set_initial_temperature(self, *args)
GridBasedPlacer_swigregister = _PcbPlacer.GridBasedPlacer_swigregister
GridBasedPlacer_swigregister(GridBasedPlacer)

# This file is compatible with both classic and new-style classes.


