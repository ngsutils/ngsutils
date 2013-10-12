import os
import subprocess

__version = None
def version():
    global __version
    if not __version:
        with open(os.path.join(os.path.realpath(os.path.dirname(__file__)), '..', 'VERSION')) as f:
            __version = f.read().strip()
            try:
                gitversion = subprocess.check_output("git show master --format='%h %ai'", cwd=os.path.join(os.path.realpath(os.path.dirname(__file__)), '..'), shell=True)
                __version = 'ngsutils-%s/%s' % (__version, gitversion.split('\n')[0].split()[0])
            except:
                pass
#    VERSION=$(echo "$GV" | awk '{print $1}')
#    echo "$(cat VERSION | sed -e 's/\n//')-$VERSION"


    return __version

