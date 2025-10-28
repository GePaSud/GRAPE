#!/bin/bash
while true; do
    CHOICE=$(dialog --clear \
        --backtitle "GRAPE Demonstration" \
        --title "Choose Example" \
        --menu "Select a simulation to run:" 15 60 4 \
        1 "Parker Solar Probe (simplified)" \
        2 "Other Example (coming soon)" \
        3 "Exit" \
        2>&1 >/dev/tty)

    clear
    case $CHOICE in
        1)
            echo "Running Parker Solar Probe example..."
            julia --project=. examples/example_ParkerSolarProbe.jl
            read -p "Press Enter to return to menu" ;;
        2)
            echo "Placeholder for future example."
            read -p "Press Enter to return to menu" ;;
        3)
            echo "Exiting GRAPE demo."
            break ;;
    esac
done
