#!/bin/bash
# Link-717 core-count test (thesis docs/handoff/l717-cores-prompt.md).
# Re-runs only the anharmonic analysis of the campaign beta-glucose ANICC+MACE4IR job from its
# checkpoint at several %NProcShared values, one job at a time, and samples the l717.exe process.
#
# usage: run_threads.sh <test_dir> <campaign_chk> <campaign_fchk> nproc [nproc ...]
set -u
T=$1; CHK=$2; FCHK=$3; shift 3
H=$(dirname "$(readlink -f "$0")")/replay_helper.py
export REPLAY_FCHK=$FCHK

sample() {  # $1 = run dir; samples l717.exe every 20 s until the job ends
    local d=$1 out=$1/samples.txt pid=""
    while [ ! -f "$d/done" ]; do
        pid=$(pgrep -n -f "l717.exe" || true)
        if [ -n "$pid" ]; then
            {
                echo "### $(date +%T) pid=$pid"
                grep -E "^(Threads|Cpus_allowed_list)" "/proc/$pid/status"
                taskset -cp "$pid"
                tr '\0' '\n' < "/proc/$pid/environ" | grep -E "^(OMP_|MKL_|GAUSS_PDEF|GAUSS_CDEF|KMP_)" || true
                echo "-- ps -L (pcpu = average over thread lifetime)"
                ps -L -o pid,tid,pcpu,stat,comm -p "$pid"
                echo "-- top -H, second frame of 2 s (instantaneous)"
                top -H -b -n 2 -d 2 -p "$pid" | awk '/^top -/{f++} f==2' | sed -n '7,40p'
            } >> "$out" 2>&1
            sleep 20
        else
            sleep 2
        fi
    done
}

for np in "$@"; do
    d=$T/np$np; mkdir -p "$d"; rm -f "$d"/done "$d"/p.chk
    while pgrep -f "l[0-9]{3}.exe" > /dev/null; do sleep 10; done  # one Gaussian job at a time
    printf '%%OldChk=%s\n%%Chk=%s\n%%mem=2GB\n%%NProcShared=%s\n#P Freq=(ReadFC,Anharm,ReadAnharm) Geom=AllCheck External="%s"\n\nDataSrc=Chk\n\n' \
        "$CHK" "$d/p.chk" "$np" "$H" > "$d/p.gjf"
    sample "$d" &
    spid=$!
    echo "$(date +%T) start np=$np" >> "$T/run_threads.log"
    (cd "$d" && REPLAY_LOG=$d/calls.log g16 < p.gjf > p.log 2>&1)
    echo "$(date +%T) end np=$np exit=$?" >> "$T/run_threads.log"
    touch "$d/done"; wait $spid
done
echo "$(date +%T) all done" >> "$T/run_threads.log"
