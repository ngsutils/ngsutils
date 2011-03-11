if [ "$1" == "" ]; then
    usage
fi

"$DIR"/build_deps.sh
if [ $? -ne 0 ]; then exit; fi

export PYTHONPATH=$PYTHONPATH:"$DIR":"$DIR"/ext

if [[ -e "$DIR"/.git && "$1" == "update" ]]; then
    cd $DIR
    if [ "$2" == "" ]; then
        branch="master"
    else
        branch="$2"
    fi
    echo "Updating from $branch branch"
    git checkout $branch
    git pull origin $branch
    "$DIR"/build_deps.sh
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
