
class Citation:
    def __init__(self, text=None, doi=None):
        # TODO: add more fields and validation as we think of them.
        self.text = text
        self.doi = doi

    # This still must be registered to be useful
    @staticmethod
    def yaml_representer(dumper, data):
        items = []
        if data.text is not None:
            items.append(('text', data.text))
        if data.doi is not None:
            items.append(('doi', data.doi))
        return dumper.represent_dict(items)
