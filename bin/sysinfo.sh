#!/bin/sh
#
# sysinfo.sh — collect hardware configuration and Java runtime information.
#
# Writes two files, both into ../results by default so that they travel with
# the other experiment results into the paper's own results folder:
#
#   results/sysinfo.txt   the full report, for the record
#   results/sysinfo.tex   the handful of facts the paper states about the
#                         machine, as LaTeX macros ready to \input
#
# Usage:
#   ./sysinfo.sh [-o OUTFILE] [-q]
#
#   -o OUTFILE   write the report here instead; the LaTeX file is then placed
#                next to it as <OUTFILE without extension>.tex
#   -q           quiet; don't echo progress to the terminal
#
# Run with sudo for complete output — DMI data (board/BIOS/RAM slots, serial
# numbers) and SMART disk details are root-only. The LaTeX file does not need
# root: everything it states comes from lscpu, /proc, lsblk and java.
#
# Plain POSIX sh, like the other scripts here, so that it runs under dash and busybox ash as well as
# bash. That rules out `-o pipefail', which run() below relied on to notice a failing command in a
# pipeline; see the comment there for what replaces it.
set -u

# Values are parsed from tool output below, so keep it in the C locale.
export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/.." 2>/dev/null && pwd) || ROOT=.

OUT=""
QUIET=0

while getopts ":o:qh" opt; do
  case "$opt" in
    o) OUT="$OPTARG" ;;
    q) QUIET=1 ;;
    h) sed -n '2,20p' "$0"; exit 0 ;;
    \?) echo "unknown option: -$OPTARG" >&2; exit 2 ;;
    :)  echo "option -$OPTARG requires an argument" >&2; exit 2 ;;
  esac
done

# A fixed name rather than a timestamped one: the paper includes the LaTeX file
# by name, and the report is regenerated together with the rest of ../results.
if [ -z "$OUT" ]; then
  mkdir -p "$ROOT/results" 2>/dev/null
  OUT="$ROOT/results/sysinfo.txt"
fi
TEXOUT="${OUT%.*}.tex"

# Scratch file for run() below, which collects a command's output before indenting it.
TMPRUN=$(mktemp "${TMPDIR:-/tmp}/sysinfo.XXXXXX") || { echo "cannot create a temporary file" >&2; exit 1; }
trap 'rm -f "$TMPRUN"' EXIT
trap 'rm -f "$TMPRUN"; exit 130' INT TERM

# ---------------------------------------------------------------- helpers ---

progress() { [ "$QUIET" -eq 1 ] || printf '  %s\n' "$*" >&2; }

section() {
  {
    printf '\n\n'
    printf '================================================================================\n'
    printf '== %s\n' "$*"
    printf '================================================================================\n'
  } >>"$OUT"
  progress "$*"
}

# run <label> <command...>
# Executes the command, appending a labelled block to the report. If the binary
# is missing the block records that instead of failing the whole run.
run() {
  run_label="$1"; shift
  {
    printf '\n--- %s\n' "$run_label"
    printf -- '$ %s\n' "$*"
  } >>"$OUT"

  if ! command -v "$1" >/dev/null 2>&1; then
    printf '    [not installed: %s]\n' "$1" >>"$OUT"
    MISSING="$MISSING $1"
    return 0
  fi

  # Capture first and indent afterwards, rather than piping the command straight into sed. A
  # pipeline reports only its *last* command's status in POSIX sh, so piped through sed a failing
  # tool would always look successful; bash's `pipefail' used to cover that, and this replaces it.
  # stderr is captured too, so tool warnings land in the report.
  if "$@" >"$TMPRUN" 2>&1; then
    run_status=0
  else
    run_status=1
  fi
  sed 's/^/    /' "$TMPRUN" >>"$OUT"
  [ "$run_status" -eq 0 ] || printf '    [command exited non-zero]\n' >>"$OUT"
}

# show <label> <file...> — dump file contents if readable
show() {
  show_label="$1"; shift
  printf '\n--- %s\n' "$show_label" >>"$OUT"
  show_found=0
  for show_f in "$@"; do
    if [ -r "$show_f" ]; then
      show_found=1
      printf '$ cat %s\n' "$show_f" >>"$OUT"
      sed 's/^/    /' "$show_f" >>"$OUT" 2>/dev/null
    fi
  done
  [ "$show_found" -eq 1 ] || printf '    [not readable: %s]\n' "$*" >>"$OUT"
}

MISSING=""

# ----------------------------------------------------------------- header ---

: >"$OUT" || { echo "cannot write $OUT" >&2; exit 1; }

{
  printf 'System & Hardware Report\n'
  printf 'Generated : %s\n' "$(date -Is 2>/dev/null || date)"
  printf 'Host      : %s\n' "$(hostname -f 2>/dev/null || hostname 2>/dev/null || echo unknown)"
  printf 'User      : %s (uid %s)\n' "${USER:-$(id -un)}" "$(id -u)"
  printf 'Kernel    : %s\n' "$(uname -srmo 2>/dev/null || uname -a)"
  if [ "$(id -u)" -ne 0 ]; then
    printf '\nNOTE: not running as root — DMI/BIOS, memory-slot and SMART data\n'
    printf '      will be incomplete. Re-run with sudo for the full picture.\n'
  fi
} >>"$OUT"

progress "writing report to $OUT"

# ------------------------------------------------------- virtualisation? ----

section "VIRTUALISATION / CONTAINER CONTEXT"
# Hardware tools report the *host* when run inside a VM or container, so
# establish this up front — it changes how the rest of the report is read.
run  "systemd-detect-virt"        systemd-detect-virt
run  "hostnamectl"                hostnamectl
show "container marker (/.dockerenv exists?)" /.dockerenv
run  "cgroup of PID 1"            cat /proc/1/cgroup

# ------------------------------------------------------- whole-machine -----

section "FULL HARDWARE INVENTORY"
run "lshw (short tree)"           lshw -short
run "lshw (detailed)"             lshw
run "hwinfo (short)"              hwinfo --short
run "inxi (full)"                 inxi -Fxxxz -c 0

section "DMI / BIOS / BASEBOARD"
run  "dmidecode (all)"            dmidecode
run  "dmidecode (memory only)"    dmidecode -t memory
show "sysfs DMI identifiers" \
     /sys/class/dmi/id/sys_vendor \
     /sys/class/dmi/id/product_name \
     /sys/class/dmi/id/product_version \
     /sys/class/dmi/id/board_vendor \
     /sys/class/dmi/id/board_name \
     /sys/class/dmi/id/bios_vendor \
     /sys/class/dmi/id/bios_version \
     /sys/class/dmi/id/bios_date \
     /sys/class/dmi/id/chassis_type

# ------------------------------------------------------------ subsystems ---

section "CPU"
run  "lscpu"                      lscpu
run  "lscpu (extended topology)"  lscpu -e
show "/proc/cpuinfo"              /proc/cpuinfo

section "MEMORY"
run  "free -h"                    free -h
show "/proc/meminfo"              /proc/meminfo
run  "swap"                       swapon --show

section "PCI DEVICES"
run "lspci -nnk (with drivers)"   lspci -nnk
run "lspci -tv (tree)"            lspci -tv

section "USB DEVICES"
run "lsusb"                       lsusb
run "lsusb -t (tree)"             lsusb -t

section "BLOCK DEVICES / STORAGE"
run "lsblk"                       lsblk -o NAME,TYPE,SIZE,ROTA,MODEL,SERIAL,TRAN,MOUNTPOINTS
run "df -hT"                      df -hT
run "mounts"                      findmnt --real

# SMART data, per physical disk.
if command -v smartctl >/dev/null 2>&1; then
  for dev in /dev/sd? /dev/nvme?n?; do
    [ -b "$dev" ] || continue
    run "smartctl -i $dev"        smartctl -i "$dev"
  done
else
  run "smartctl"                  smartctl --version
fi

section "NETWORK"
run "interfaces (ip -br)"         ip -br address
run "links"                       ip -br link
run "routes"                      ip route
if command -v ethtool >/dev/null 2>&1; then
  for nic in $(ls /sys/class/net 2>/dev/null); do
    [ "$nic" = lo ] && continue
    run "ethtool $nic"            ethtool "$nic"
    run "ethtool -i $nic"         ethtool -i "$nic"
  done
else
  run "ethtool"                   ethtool --version
fi

section "GRAPHICS"
run "glxinfo (OpenGL)"            glxinfo -B
run "nvidia-smi"                  nvidia-smi
show "DRM cards"                  /sys/class/drm/card0/device/uevent

section "SENSORS / POWER"
run "sensors"                     sensors
run "upower"                      upower -d

section "OPERATING SYSTEM"
show "/etc/os-release"            /etc/os-release
run  "uname -a"                   uname -a
run  "uptime"                     uptime
run  "loaded modules"             lsmod

# ------------------------------------------------------------------ Java ---

section "JAVA VIRTUAL MACHINE"

printf '\n--- environment\n' >>"$OUT"
printf '    JAVA_HOME=%s\n'  "${JAVA_HOME:-<unset>}"  >>"$OUT"
printf '    JDK_HOME=%s\n'   "${JDK_HOME:-<unset>}"   >>"$OUT"
printf '    JAVA_OPTS=%s\n'  "${JAVA_OPTS:-<unset>}"  >>"$OUT"
printf '    java on PATH: %s\n' "$(command -v java || echo '<none>')" >>"$OUT"

# Version of the default runtime and compiler. Note: `java -version` writes to
# stderr, which run() already folds into the report.
run "java -version"               java -version
run "java --version"              java --version
run "javac -version"              javac -version

# The authoritative dump of what this VM actually is: vendor, VM name and
# variant (HotSpot/OpenJ9), heap sizing, boot classpath, os.arch.
run "java -XshowSettings (all properties)" \
    java -XshowSettings:all -version

# Only the flags that differ from their defaults — ergonomics (GC choice, heap
# sizing, thread counts) derived by the VM from the hardware it found. The full
# list is ~1500 lines; drop the grep below if you want all of it.
run "java VM flags (non-default / ergonomic)" \
    sh -c "java -XX:+PrintFlagsFinal -version 2>/dev/null | grep -Ev '\{default\}' "
run "available processors / heap seen by the VM" \
    java -XshowSettings:system -version

# All JVMs installed on the machine, not just the one first on PATH.
printf '\n--- installed JVMs\n' >>"$OUT"
for d in /usr/lib/jvm /usr/java /opt/java /Library/Java/JavaVirtualMachines "$HOME/.sdkman/candidates/java"; do
  [ -d "$d" ] || continue
  printf '$ ls -1 %s\n' "$d" >>"$OUT"
  ls -1 "$d" 2>/dev/null | sed 's/^/    /' >>"$OUT"
done

# Report the version of every JVM found, not only the default one.
printf '\n--- version of each installed JVM\n' >>"$OUT"
found_jvm=0
for jbin in /usr/lib/jvm/*/bin/java /usr/java/*/bin/java /opt/java/*/bin/java \
            "$HOME"/.sdkman/candidates/java/*/bin/java; do
  [ -x "$jbin" ] || continue
  found_jvm=1
  printf '$ %s -version\n' "$jbin" >>"$OUT"
  "$jbin" -version 2>&1 | sed 's/^/    /' >>"$OUT"
done
[ "$found_jvm" -eq 1 ] || printf '    [no JVMs found in the usual locations]\n' >>"$OUT"

run "update-alternatives (java)"  update-alternatives --display java
run "sdkman current"              sdk current java

# Running JVM processes and their VM identification.
run "running JVMs (jps)"          jps -lvm
if command -v jcmd >/dev/null 2>&1; then
  printf '\n--- jcmd VM.version for running JVMs\n' >>"$OUT"
  # Column 1 of `jcmd -l` is the PID; skip jcmd's own entry.
  jcmd -l 2>/dev/null | while read -r pid rest; do
    case "$rest" in *jdk.jcmd*) continue ;; esac
    printf '$ jcmd %s VM.version   (%s)\n' "$pid" "$rest" >>"$OUT"
    jcmd "$pid" VM.version 2>&1 | sed 's/^/    /' >>"$OUT"
  done
fi

# Package-manager view, useful for auditing which JDK packages are installed.
run "dpkg java packages"          sh -c "dpkg -l | grep -Ei 'jdk|jre|java'"
run "rpm java packages"           sh -c "rpm -qa | grep -Ei 'jdk|jre|java'"

# --------------------------------------------------------- LaTeX extract ---
#
# The report above is for the record; what the paper actually states about the
# machine is a handful of numbers. They are collected here into LaTeX macros so
# that no one has to copy them by hand -- see the "Performance" section of the
# Genestrip-FT paper, which mirrors the hardware description of the first one.
#
# Values are queried from the system directly rather than parsed back out of
# the report, which keeps this independent of the report's formatting. Anything
# that cannot be determined becomes \sysUnknown, which typesets as a bold "??"
# and is therefore impossible to miss in the PDF.

section "LATEX EXTRACT FOR THE PAPER"

# Escapes the characters that would otherwise break or silently alter LaTeX.
tex_escape() {
  sed -e 's/\\/\\textbackslash{}/g' \
      -e 's/\([&%$#_{}]\)/\\\1/g' \
      -e 's/~/\\textasciitilde{}/g' \
      -e 's/\^/\\textasciicircum{}/g'
}

lscpu_field() {
  lscpu 2>/dev/null | sed -n "s/^$1: *//p" | head -1 | sed 's/ *$//'
}

# Tools do not always leave an unknown field empty -- lscpu prints "-" for the
# CPU model inside a VM, dmidecode has its own set of placeholders. Treating
# those as values would put a literal "-" into the paper.
sane() {
  case "$(printf '%s' "${1:-}" | sed 's/^ *//; s/ *$//')" in
    ''|-|unknown|Unknown|n/a|N/A|None|'Not Specified'|'To Be Filled By O.E.M.') printf '' ;;
    *) printf '%s' "$1" ;;
  esac
}

# --- CPU
cpu_model=$(sane "$(lscpu_field 'Model name')")
[ -n "$cpu_model" ] || cpu_model=$(sed -n 's/^model name[[:space:]]*: *//p' /proc/cpuinfo 2>/dev/null | head -1)
# "Intel(R) Xeon(R) W-2255 CPU @ 3.70GHz" reads better without the trademarks.
cpu_model=$(printf '%s' "$cpu_model" | sed -e 's/(R)//g' -e 's/(TM)//g' -e 's/  */ /g' -e 's/ *$//')

# The nominal clock, taken from the model name where the vendor put it there
# and from lscpu otherwise. lscpu reports MHz, hence the division.
cpu_ghz=$(printf '%s' "$cpu_model" | sed -n 's/.*@ *\([0-9.]*\) *GHz.*/\1/p')
if [ -z "$cpu_ghz" ]; then
  cpu_mhz=$(lscpu_field 'CPU max MHz')
  [ -n "$cpu_mhz" ] || cpu_mhz=$(lscpu_field 'CPU MHz')
  [ -n "$cpu_mhz" ] && cpu_ghz=$(awk -v m="$cpu_mhz" 'BEGIN{printf "%.2f", m/1000}')
fi

# Physical cores, which is what a paper quotes -- not the logical CPUs that
# hyper-threading multiplies them into. Three ways to arrive at the number,
# because lscpu leaves different fields blank depending on the platform.
sockets=$(sane "$(lscpu_field 'Socket(s)')")
per_socket=$(sane "$(lscpu_field 'Core(s) per socket')")
threads=$(sane "$(lscpu_field 'CPU(s)')")
per_core=$(sane "$(lscpu_field 'Thread(s) per core')")
cores=""
if [ -n "$sockets" ] && [ -n "$per_socket" ]; then
  cores=$((sockets * per_socket))
elif [ -n "$threads" ] && [ -n "$per_core" ] && [ "$per_core" -gt 0 ] 2>/dev/null; then
  cores=$((threads / per_core))
else
  # Last resort: count the distinct (socket, core) pairs the kernel reports.
  cores=$(awk -F: '/^physical id/{p=$2} /^core id/{print p":"$2}' /proc/cpuinfo 2>/dev/null \
          | sort -u | grep -c .)
  [ "${cores:-0}" -gt 0 ] 2>/dev/null || cores=""
fi

# The worker threads the Genestrip goals actually run with: pom.xml sets
# gs.threads to -1, which Genestrip reads as one thread per available processor
# less one. Deriving it here keeps the number the paper quotes tied to the
# machine that produced the measurements.
workers=""
if [ -n "$threads" ] && [ "$threads" -gt 1 ] 2>/dev/null; then
  workers=$((threads - 1))
fi

# --- memory
mem_kb=$(sed -n 's/^MemTotal: *\([0-9]*\).*/\1/p' /proc/meminfo 2>/dev/null)
ram_gib=""
ram_nominal=""
if [ -n "$mem_kb" ]; then
  # MemTotal excludes what the firmware and the kernel reserve, so a 128 GB
  # machine reports about 125 GiB. Rounding up to the next multiple of 8
  # recovers the size one would quote in a paper; \sysRamGiB keeps the raw
  # figure so the rounding can be checked.
  ram_gib=$(( (mem_kb + 1048575) / 1048576 ))
  ram_nominal=$(( ((ram_gib + 7) / 8) * 8 ))
fi

swap_kb=$(sed -n 's/^SwapTotal: *\([0-9]*\).*/\1/p' /proc/meminfo 2>/dev/null)
if [ -z "$swap_kb" ]; then
  swap_state=""
elif [ "$swap_kb" -eq 0 ]; then
  swap_state="disabled"
else
  swap_state="enabled ($(( (swap_kb + 1048575) / 1048576 )) GiB)"
fi

# --- disk holding the experiment data
disk_size=""
disk_type=""
disk_model=""
data_dir="$ROOT/data"
[ -d "$data_dir" ] || data_dir="$ROOT"
# findmnt resolves a bind mount to its underlying source, which df does not.
part=$(findmnt -no SOURCE --target "$data_dir" 2>/dev/null | head -1)
[ -n "${part:-}" ] || part=$(df -P "$data_dir" 2>/dev/null | awk 'NR==2{print $1}')
# A bind mount is reported as "/dev/sda2[/sub/dir]"; only the device counts.
part=${part%%[*}
if [ -n "${part:-}" ] && command -v lsblk >/dev/null 2>&1; then
  # The partition's parent is the physical disk; on a whole-device mount there
  # is no parent and the device itself is the disk.
  dev=$(lsblk -no PKNAME "$part" 2>/dev/null | head -1)
  [ -n "$dev" ] || dev=$(lsblk -dno KNAME "$part" 2>/dev/null | head -1)
  if [ -n "$dev" ] && [ -b "/dev/$dev" ]; then
    bytes=$(lsblk -dbno SIZE "/dev/$dev" 2>/dev/null | head -1)
    rota=$(lsblk -dno ROTA "/dev/$dev" 2>/dev/null | head -1 | tr -d ' ')
    tran=$(lsblk -dno TRAN "/dev/$dev" 2>/dev/null | head -1 | tr -d ' ')
    disk_model=$(lsblk -dno MODEL "/dev/$dev" 2>/dev/null | head -1 | sed 's/ *$//')
    if [ -n "${bytes:-}" ]; then
      # Manufacturers count in powers of ten, so that is how a disk is quoted.
      disk_size=$(awk -v b="$bytes" 'BEGIN{
        if (b >= 1e12) printf "%.3g TB", b/1e12; else printf "%.0f GB", b/1e9 }')
    fi
    case "$rota" in
      0) disk_type="solid state" ;;
      1) disk_type="rotational" ;;
    esac
    [ "$tran" = nvme ] && [ "$disk_type" = "solid state" ] && disk_type="NVMe solid state"
  fi
fi

# --- operating system
os_name=$(sed -n 's/^PRETTY_NAME="\(.*\)"$/\1/p' /etc/os-release 2>/dev/null | head -1)
[ -n "$os_name" ] || os_name=$(uname -sr 2>/dev/null)
kernel=$(uname -r 2>/dev/null)

# --- Java: the runtime that actually executes the goals, i.e. the one on PATH
java_props=$(java -XshowSettings:properties -version 2>&1)
java_prop() { printf '%s' "$java_props" | sed -n "s/^ *$1 = *//p" | head -1 | sed 's/ *$//'; }
java_vendor=$(java_prop 'java.vendor')
java_version=$(java_prop 'java.version')
java_vm=$(java_prop 'java.vm.name')
java_vm_version=$(java_prop 'java.vm.version')

# --- write the file
# def <macro> <value> -- emits \newcommand, falling back to \sysUnknown
def() {
  if [ -n "${2:-}" ]; then
    printf '\\newcommand{\\%s}{%s}\n' "$1" "$(printf '%s' "$2" | tex_escape)" >>"$TEXOUT"
  else
    printf '\\newcommand{\\%s}{\\sysUnknown}\n' "$1" >>"$TEXOUT"
    MISSING_TEX="$MISSING_TEX $1"
  fi
}

MISSING_TEX=""
: >"$TEXOUT" || { echo "cannot write $TEXOUT" >&2; exit 1; }
{
  printf '%% Hardware and runtime of the machine the experiments were executed on.\n'
  printf '%% Generated by ft-db-exp2/bin/sysinfo.sh on %s -- do not edit by hand.\n' \
         "$(date -Is 2>/dev/null || date)"
  printf '%% Host: %s\n' "$(hostname -f 2>/dev/null || hostname 2>/dev/null || echo unknown)"
  printf '%% Copy this file to wherever the other results are read from -- it belongs\n'
  printf '%% with them; a LaTeX consumer inputs it and falls back to bold "??" markers\n'
  printf '%% for anything it does not find.\n'
  printf '\\newcommand{\\sysUnknown}{\\textbf{??}}\n'
} >>"$TEXOUT"

def sysCpuModel     "$cpu_model"
def sysCpuGHz       "${cpu_ghz:-}"
def sysCores        "${cores:-}"
def sysThreads      "${threads:-}"
def sysWorkers      "${workers:-}"
# The threads actually doing work: the workers plus the one thread that reads and parses the input.
# Emitted as a finished number rather than left to LaTeX arithmetic in the paper. The earlier
# version passed a \newif flag instead, so that the paper could tell whether \sysWorkers was a value
# it could compute with -- and that broke the document: TeX counts the \if... token of a
# `\newif\ifX' while it skips the branch that \newif sits in, once \ifX is defined, so the guard
# `\ifdefined\ifX\else\newif\ifX\fi' swallowed its own \fi and the run ended in "Incomplete
# \ifdefined". Nothing here needs a flag; an unknown value simply becomes \sysUnknown like any other.
total_threads=""
if [ -n "${workers:-}" ]; then
  total_threads=$((workers + 1))
fi
def sysTotalThreads "${total_threads:-}"
def sysRamGB        "${ram_nominal:-}"
def sysRamGiB       "${ram_gib:-}"
def sysSwap         "${swap_state:-}"
def sysDiskSize     "${disk_size:-}"
def sysDiskType     "${disk_type:-}"
def sysDiskModel    "${disk_model:-}"
def sysOs           "${os_name:-}"
def sysKernel       "${kernel:-}"
def sysJavaVendor   "${java_vendor:-}"
def sysJavaVersion  "${java_version:-}"
def sysJavaVm       "${java_vm:-}"
def sysJavaVmVersion "${java_vm_version:-}"

{
  printf '\n--- %s\n' "$TEXOUT"
  sed 's/^/    /' "$TEXOUT"
  if [ -n "$MISSING_TEX" ]; then
    printf '\n    NOT DETERMINED:%s\n' "$MISSING_TEX"
    printf '    These appear as bold "??" in the PDF and must be filled in by hand.\n'
  fi
} >>"$OUT"

if [ "$QUIET" -ne 1 ]; then
  [ -n "$MISSING_TEX" ] && printf '  WARNING: not determined:%s\n' "$MISSING_TEX" >&2
fi

# ------------------------------------------------------------------ tail ---

section "REPORT SUMMARY"
if [ -n "$MISSING" ]; then
  printf '\nTools not installed (sections above are incomplete):\n' >>"$OUT"
  printf '%s\n' "$MISSING" | tr ' ' '\n' | grep -v '^$' | sort -u | sed 's/^/    /' >>"$OUT"
  printf '\nOn Debian/Ubuntu install them with:\n' >>"$OUT"
  printf '    sudo apt-get install -y lshw dmidecode pciutils usbutils \\\n' >>"$OUT"
  printf '        smartmontools ethtool lm-sensors inxi hwinfo mesa-utils\n' >>"$OUT"
else
  printf '\nAll queried tools were present.\n' >>"$OUT"
fi
printf '\nEnd of report.\n' >>"$OUT"

if [ "$QUIET" -ne 1 ]; then
  printf '\nReport written to : %s (%s)\n' "$OUT" "$(du -h "$OUT" | cut -f1)" >&2
  printf 'LaTeX macros      : %s\n' "$TEXOUT" >&2
fi
