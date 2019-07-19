import collections
import re

class DRUFile():
    """Design rule options for Eagle.  This class loads configuration values from
    a :code:`dru` file and exposes them via :meth:`get_value()` and as
    attributes on the object.

    It converts all values with units provided in the file into millimeters.
    """

    def __init__(self, stream, filename=None):
        """
        Open a DRU file and load its contents.
        
        :param stream: File-like object to load.
        """

        self.values = collections.OrderedDict()
        self.open(stream, filename=filename);
        
    def open(self, stream, filename=None):
        """
        Open a DRU file and load its contents.
        
        :param stream: File-like object to load.
        """

        for l in stream.readlines():
            m = re.match("^(\w+)(\[(\w+)\])? = (.*)$", l);
            assert m is not None, "Unexpected line format in '{}': {}".format(filename,l)
            key = m.group(1)
            dictkey = m.group(3)
            value = m.group(4)

            values = value.split(" ")
            if key in ["description", "layerSetup"]:
                values = [value]
            else:
                values = value.split(" ")
                for i in range(0, len(values)):
                    m = re.match("^(-?\d+(\.\d+)?)(\w+)?$", values[i]);
                    assert m is not None, "Unknown value format in '{}': {}".format(filename,values[i])
                    n = m.group(1)
                    units = m.group(3)
                    if units is None:
                        values[i] = float(n)
                    elif units == "mm":
                        values[i] = float(n)
                    elif units == "mil":
                        values[i] = float(n)*0.0254
                    else:
                        assert False, "Unknown unit in '{}': {}".format(filename,units)

            m = re.match("(\w+)(\[(\w+)\])?", key);
            assert m is not None, "Unknown key format in '{}': {}".format(filename,key)
            key = m.group(1);

            if len(values) == 1:
                value = values[0]
            else:
                value = values

            if dictkey is not None:
                self.values.setdefault(key, {})[dictkey] = value
            else:
                self.values[key] = value

        for k in self.values:
            setattr(self, k, self.values[k])



    def get_value(self, key):
        """
        Get a value by its key.
        
        :param key: The key.
        :returns: The corresponding value.  It may be a single, array, or map, depending on the contents of the DRU file.
        :rtype: varies.
        """
        return self.values[key]
