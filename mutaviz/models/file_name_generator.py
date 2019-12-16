from random import choice
import string
import os


class FileNameGenerator:
    def random(self, extension=None, path=os.getcwd()):
        file_path_name = (path + "/") + "".join(choice(string.ascii_letters) for x in range(20)) + ".%s" % extension
        return os.path.abspath(file_path_name)
