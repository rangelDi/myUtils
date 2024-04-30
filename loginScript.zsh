#!/bin/zsh
chrono -p
echo "Enter 'y' to continue..."
while true; do
    # Prompt the user for input
    read response

    # Check if the response is 'y' or 'Y'
    if [[ "$response" == "y" || "$response" == "Y" ]]; then
        break  # Exit the loop if the response is 'y' or 'Y'
    else
        echo "Enter 'y' to continue..."
    fi
done
