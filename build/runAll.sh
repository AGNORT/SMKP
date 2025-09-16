#!/bin/bash

# path to the executable file
EXEC="./cbp"

# path to the instances files
INSTANCE_DIR=""

# path to the results files
RESULT_DIR=""
mkdir -p "$RESULT_DIR"

# run all instances
for instance in "$INSTANCE_DIR"/*.txt; do
    filename=$(basename "$instance" .txt)
    
    result="$RESULT_DIR/$filename.txt"
    
    echo "Running: $EXEC $instance $result"
    $EXEC "$instance" "$result"
done

