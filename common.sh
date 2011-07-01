if [ "$1" == "" ]; then
    usage
fi

REAL=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
DIR=`dirname "$REAL"`

# "$DIR"/build_deps.sh
# if [ $? -ne 0 ]; then exit; fi

export PYTHONPATH=$PYTHONPATH:"$DIR":"$DIR"/ext

if [[ -e "$DIR"/.git && "$1" == "update" ]]; then
    cd "$DIR"
    
    OLD_PYSAM_VER=`grep 'PYSAM_VER' build_deps.sh`
    
    if [ "$2" != "" ]; then
        echo "Updating from $2 branch"
        git checkout $2
        git pull origin $2
    else
        echo "Updating from current branch"
        git pull
    fi

    NEW_PYSAM_VER=`grep 'PYSAM_VER' build_deps.sh`

    if [ "$OLD_PYSAM_VER" != "$NEW_PYSAM_VER" ]; then
        rm -rf ext/pysam
        echo "Updating dependencies..."
        ./build_deps.sh
    fi

    exit 0
fi


if [ "$1" == "help" ]; then
    if [ "$2" == "" ]; then
        usage
    fi
    
    action=$PREFIX$2.py
    
    if [ ! -e "$DIR"/$SUBDIR/$action ]; then
        action=$PREFIX$2.sh
        if [ ! -e "$DIR"/$SUBDIR/$action ]; then
            echo "Unknown command '$2'"
            exit 1
        fi
    fi
    "$DIR"/$SUBDIR/$action -h
else
    action=$PREFIX$1.py
    
    if [ ! -e "$DIR"/$SUBDIR/$action ]; then
        action=$PREFIX$1.sh
        if [ ! -e "$DIR"/$SUBDIR/$action ]; then
            echo "Unknown command '$1'"
            exit 1
        fi
    fi
    shift

    ARGS=()
    i=0
    for arg in "$@"; do
        ARGS[$i]="$arg"
        ((++i))
    done
    
    exec "$DIR"/$SUBDIR/$action "${ARGS[@]}"
fi
