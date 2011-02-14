if [ "$1" == "" ]; then
    usage
fi

"$DIR"/build_deps.sh
if [ $? -ne 0 ]; then exit; fi

export PYTHONPATH=$PYTHONPATH:"$DIR":"$DIR"/ext

if [ "$1" == "update" ]; then
    cd $DIR
    git pull
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

    "$DIR"/$SUBDIR/$action $@
fi
