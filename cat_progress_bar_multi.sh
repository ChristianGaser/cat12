#!/usr/bin/env bash
#
# PURPOSE: Displays progress bars (multi-job or single-job).
#          Adopted from T1Prep (scripts/progress_bar_multi.sh).
#
# USAGE:
#   # Single job (no progress dir):
#   cat_progress_bar_multi.sh 1 "" CURRENT TOTAL [Label Width FailedFlag]
#
#   # Multi job:
#   cat_progress_bar_multi.sh n_jobs PROGRESS_DIR [Width Label] --multi
#
#   # Multi job (overall bar only):
#   cat_progress_bar_multi.sh n_jobs PROGRESS_DIR [Width Label]
#
# ______________________________________________________________________
#
# Christian Gaser
# Structural Brain Mapping Group (https://neuro-jena.github.io)
# Departments of Neurology and Psychiatry
# Jena University Hospital
# ______________________________________________________________________
# $Id$

# ----------------------------------------------------------------------
# global parameters
# ----------------------------------------------------------------------
n_jobs=$1
PROGRESS_DIR=$2
width=$3
job_name=$4
overall=1
[ "$5" = "--multi" ] && overall=0

# Colors & cursor fallback
if [ -t 1 ] && [ -n "$TERM" ] && command -v tput >/dev/null 2>&1 && \
   tput sgr0 >/dev/null 2>&1; then
  BOLD=$(tput bold)
  NC=$(tput sgr0) # Reset all attributes
  DEFAULT_COLOR=$(tput setaf 2)   # green
  FAIL_COLOR=$(tput setaf 1)      # red
  move_up() { tput cuu "$1"; }
elif [ -t 1 ]; then
  # terminal without usable terminfo entry (e.g. TERM=dumb): use plain ANSI
  # codes and skip the colors
  move_up() { printf '\033[%dA' "$1"; }
else
  move_up() { :; }
fi

REFRESH_INTERVAL=1.0
JOB_ID="jobrun_$(date +%s%N)"

# ----------------------------------------------------------------------
# Single-job mode
# ----------------------------------------------------------------------
single_job_bar() {
  local current=$1 total=$2 label=$3 width=$4 failed=$5
  [ -z "$label" ] && label="Progress"
  [ -z "$width" ] && width=40
  [ -z "$failed" ] && failed=0

  # Track elapsed time across separate script invocations to estimate ETA.
  local state_dir="/tmp/cat_progress_single"
  mkdir -p "$state_dir"
  local state_file="$state_dir/${PPID}.state"
  local now start prev_total prev_current
  now=$(date +%s)

  if [[ "$current" -le 1 || ! -f "$state_file" ]]; then
    start=$now
    printf "%s,%s,%s\n" "$start" "$total" "$current" > "$state_file"
  else
    IFS=',' read -r start prev_total prev_current < "$state_file"
    [ -z "$start" ] && start=$now
    [ -z "$prev_total" ] && prev_total=$total
    [ -z "$prev_current" ] && prev_current=0
    # Reset timing when a new run starts or total changed.
    if [[ "$current" -le "$prev_current" || "$total" -ne "$prev_total" ]]; then
      start=$now
    fi
    printf "%s,%s,%s\n" "$start" "$total" "$current" > "$state_file"
  fi

  local percent=$(( current * 100 / total ))
  local filled=$(( width * current / total ))
  local empty=$(( width - filled ))

  local COLOR=$DEFAULT_COLOR
  [ "$failed" -ne 0 ] && COLOR=$FAIL_COLOR

  local eta=""
  if [[ "$current" -gt 1 && "$current" -lt "$total" ]]; then
    local elapsed remaining
    elapsed=$((now - start))
    if [[ "$elapsed" -gt 0 ]]; then
      remaining=$((elapsed * (total - current) / (current - 1)))
      eta=$(printf "ETA %02d:%02d:%02d" \
        $((remaining / 3600)) $(((remaining / 60) % 60)) $((remaining % 60)))
    fi
  fi

  local bar=$(printf "%${filled}s" "" | awk '{gsub(/ /,"█"); print}')
  bar+=$(printf "%${empty}s")

  if [ -t 1 ]; then
    printf "|${COLOR}%-${width}s${NC}| %3d%% %s %s\r" "$bar" "$percent" "$label" "$eta"
    if [ "$current" -eq "$total" ]; then
      rm -f "$state_file"
      echo ""  # newline at end
    fi
  else
    # Non-TTY: only final output
    if [ "$current" -eq "$total" ]; then
      rm -f "$state_file"
      echo "|${COLOR}${bar}${NC}| ${percent}% $label"
    fi
  fi

}

# ----------------------------------------------------------------------
# Multi-job mode
# ----------------------------------------------------------------------
multi_job_bar() {
  [ -z "$PROGRESS_DIR" ] && PROGRESS_DIR="/tmp/progress_bars/$JOB_ID"
  [ -z "$width" ] && width=40
  [ -z "$job_name" ] && job_name="Job"

  local start_time=$(date +%s)

  if [ -t 1 ]; then
    if [ $overall -eq 0 ]; then
      for ((i=0; i<n_jobs; i++)); do
        echo ""
      done
    else
      echo ""
    fi
  fi

  while true; do
    if [ -t 1 ]; then
      if [ $overall -eq 0 ]; then
        move_up "$n_jobs"
      else
        move_up 1
      fi
    fi

    all_done=true

    if [ $overall -eq 1 ]; then
      total_done=0
      total_items_sum=0
    fi

    for ((i=0; i<n_jobs; i++)); do
      file="$PROGRESS_DIR/job${i}.progress"
      if [[ -f "$file" ]]; then
        progress=$(cat "$file")
        done_items=$(echo "$progress" | cut -d/ -f1)
        total_items=$(echo "$progress" | cut -d/ -f2)

        # Initialise missing or malformed files safely
        [ -z "$done_items" ] && done_items=0
        [ -z "$total_items" ] && total_items=1

        if [ $overall -eq 1 ]; then
          total_done=$((total_done + done_items))
          total_items_sum=$((total_items_sum + total_items))
        else
          percent=$((100 * done_items / total_items))
          filled=$((width * done_items / total_items))
          unfilled=$((width - filled))
          bar=$(printf "%${filled}s" "" | awk '{gsub(/ /,"█"); print}')
          bar+=$(printf "%${unfilled}s")
          jobnumber=$((i+1))

          if [ $n_jobs -gt 1 ]; then
            printf -v jobstr "%2d" "$jobnumber"
          fi
          printf -v percent_str "%3d" "$percent"

          status_file="$PROGRESS_DIR/job${i}.status"
          status=0
          [[ -f "$status_file" ]] && status=$(cat "$status_file")
          if [[ "$status" -ne 0 ]]; then
            BAR_COLOR=$FAIL_COLOR
            REFRESH_INTERVAL=0.1
          else
            BAR_COLOR=$DEFAULT_COLOR
            REFRESH_INTERVAL=2.0
          fi

          if [ -t 1 ]; then
            echo -ne "|${BAR_COLOR}${bar}${NC}| ${percent_str}% ${job_name} ${jobstr}\n"
          else
            if [[ "$done_items" -eq "$total_items" ]]; then
              echo "|${BAR_COLOR}${bar}${NC}| ${percent_str}% ${job_name} ${jobstr}"
            fi
          fi
        fi

        [[ "$done_items" -lt "$total_items" ]] && all_done=false
    else
      # Initialise progress file with 0/TOTAL = 1 to avoid div-by-0
      echo "0/1" > "$PROGRESS_DIR/job${i}.progress"
      if [ -t 1 ] && [ $overall -eq 0 ]; then
        echo -ne "Job $((i+1)): [initialising...]\n"
      fi
      all_done=false
      fi
    done

    if [ $overall -eq 1 ]; then
      if [ "$total_items_sum" -gt 0 ]; then
        percent=$((100 * total_done / total_items_sum))
        filled=$((width * total_done / total_items_sum))

        # Compute ETA
        if [ $total_done -gt 0 ] && [ $total_done -lt $total_items_sum ]; then
          now=$(date +%s)
          elapsed=$((now - start_time))
          # Remaining time = elapsed * (total/so_far - 1)
          remaining=$((elapsed * (total_items_sum - total_done) / total_done))

          # Format as HH:MM:SS
          ETA=$(printf "%02d:%02d:%02d" $((remaining/3600)) $(( (remaining/60)%60 )) $((remaining%60)))
          ETA_str="ETA"
        else
          ETA=""
          ETA_str=""
        fi
      else
        percent=0
        filled=0
        ETA=""
        ETA_str=""
      fi

      unfilled=$((width - filled))
      bar=$(printf "%${filled}s" "" | awk '{gsub(/ /,"█"); print}')
      bar+=$(printf "%${unfilled}s")
      printf "|${DEFAULT_COLOR}${bar}${NC}| %4d%% %s %s %s\n" "$percent" "$job_name" "$ETA_str" "$ETA"
    fi

    $all_done && break
    sleep "$REFRESH_INTERVAL"
  done
}

# ----------------------------------------------------------------------
# Main logic
# ----------------------------------------------------------------------
main() {
  if [ -z "$n_jobs" ]; then
    help
    exit 1
  fi

  # Hide cursor
  if [ -t 1 ]; then
    echo -ne "\033[?25l"
    trap 'echo -ne "\033[?25h"' EXIT INT TERM  # Ensure cursor restored on exit
  fi

  if [ "$n_jobs" -eq 1 ] && [ -z "$PROGRESS_DIR" ]; then
    # Single-job mode
    single_job_bar "$3" "$4" "$5" "$6" "$7"
  else
    # Multi-job mode
    multi_job_bar
  fi
}

# ----------------------------------------------------------------------
# Help
# ----------------------------------------------------------------------
help() {
cat <<__EOM__

${BOLD:-}USAGE:${NC}
  Single job:
    $(basename "$0") 1 "" CURRENT TOTAL [Label Width FailedFlag]

  Multi job:
    $(basename "$0") n_jobs PROGRESS_DIR [Width Label] --multi

  Multi job (overall bar only):
    $(basename "$0") n_jobs PROGRESS_DIR [Width Label]

${BOLD:-}OPTIONS:${NC}
  CURRENT      Current progress value (single job).
  TOTAL        Maximum progress value (single job).
  Label        Optional name of process (default "Progress"/"Job").
  Width        Width of progress bar (default 40).
  FailedFlag   Set to 1 to render the bar in red (single job only).

  n_jobs       Number of jobs to monitor in parallel (multi job).
  PROGRESS_DIR Directory containing job*.progress and job*.status files.
  --multi      Show (multiple) progress bars for each job.

${BOLD:-}EXAMPLES:${NC}

  # --- Single job ---
  # Progress 5 of 10 (50%) for a task called "SkullStrip":
  $(basename "$0") 1 "" 5 10 "SkullStrip"

  # --- Multi job ---
  # 3 parallel tasks, each job writes its own progress file:
  n_jobs=3
  PROGRESS_DIR="/tmp/progress_bars"
  mkdir -p "\$PROGRESS_DIR"
  rm -f "\$PROGRESS_DIR"/*.progress

  # Dummy command
  COMMAND="sleep"  # Replace with your actual command

  # Arguments per job (simulate multiple tasks per job)
  ARG_LIST=(
    "1 2 3 4 5"
    "a b c"
    "x y z w"
  )

  # Start progress monitor
  $(dirname "\$0")/$(basename "$0") "\$n_jobs" "\$PROGRESS_DIR" 50 "Stage" &
  progress_pid=\$!

  # Launch jobs
  for ((i=0; i<n_jobs; i++)); do
    (
      args=(\${ARG_LIST[\$i]})
      total=\${#args[@]}
      for ((k=0; k<total; k++)); do
        echo "\$((k+1))/\$total" > "\$PROGRESS_DIR/job\${i}.progress"
        \$COMMAND 0.5  # Simulate work
      done
    ) &
  done

  wait  # Wait for all jobs

  echo -e "\\n✅ All jobs completed!"
  rm -r "\$PROGRESS_DIR"

${BOLD:-}Author:${NC}
  Christian Gaser (christian.gaser@uni-jena.de)

__EOM__
}

# ----------------------------------------------------------------------
# call main program
# ----------------------------------------------------------------------

main ${1+"$@"}
