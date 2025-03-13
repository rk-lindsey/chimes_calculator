#!/bin/bash

# Process all cluster types (2b, 3b, 4b)
for body in 2b 3b 4b; do
    echo "Processing ${body} clusters..."
    
    # Find all cluster files for this body type
    find . -maxdepth 1 -name "*.[0-9]*.${body}_clusters.txt" -print0 | while IFS= read -r -d $'\0' file; do
        # Extract components from filename
        filename=$(basename "$file")
        frame=$(echo "$filename" | cut -d. -f1)
        rank=$(echo "$filename" | cut -d. -f2)
        
        # Store in temporary frame list
        echo "$rank $file" >> "tmp_${body}_${frame}.list"
    done

    # Process each frame's temporary file
    for framefile in tmp_${body}_*.list; do
        if [ -f "$framefile" ]; then
            frame=${framefile#tmp_${body}_}
            frame=${frame%.list}
            
            echo "  Combining frame $frame for ${body}"
            
            # Sort files by rank and concatenate
            sort -n -k1 "$framefile" | awk '{print $2}' | xargs cat > "${frame}.${body}_combined.txt"
            
            # Remove original files and corresponding lists if concatenation succeeded
            if [ $? -eq 0 ]; then
                echo "  Removing original files:"
                # Generate list of files to remove (both clusters and lists)
                awk '{
                    cluster_file = $2; 
                    list_file = cluster_file; 
                    sub(/_clusters\.txt$/, "_list.txt", list_file); 
                    print cluster_file; 
                    print list_file
                }' "$framefile" | xargs rm -v
            else
                echo "  Error concatenating files, preserving originals"
            fi
            
            # Clean up temporary file
            rm "$framefile"
        fi
    done
done

echo -e "\nCombination complete. Created files:"
ls -l *b_combined.txt