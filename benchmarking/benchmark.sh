#!/bin/bash

# C++ Programm Benchmark Script für M1 Mac mit Xcode Tools
# Verwendung: ./benchmark.sh ./<programm_name> [anzahl_läufe]

set -e

# Farben für Output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Funktionen
print_usage() {
    echo -e "${BLUE}Verwendung:${NC} $0 <programm_pfad> [anzahl_läufe]"
    echo -e "${BLUE}Beispiel:${NC} $0 ./mein_programm 5"
    echo ""
    echo -e "${BLUE}Optionen:${NC}"
    echo "  programm_pfad    Pfad zum kompilierten C++ Programm"
    echo "  anzahl_läufe     Anzahl der Benchmark-Läufe (Standard: 3)"
}

print_header() {
    echo -e "${GREEN}================================${NC}"
    echo -e "${GREEN} C++ Programm Benchmark (M1 Mac)${NC}"
    echo -e "${GREEN}================================${NC}"
    echo ""
}

check_dependencies() {
    echo -e "${YELLOW}Überprüfe Abhängigkeiten...${NC}"
    
    # Prüfe ob Xcode Command Line Tools installiert sind
    if ! command -v xcrun &> /dev/null; then
        echo -e "${RED}Fehler: Xcode Command Line Tools nicht gefunden!${NC}"
        echo "Installiere sie mit: xcode-select --install"
        exit 1
    fi
    
    # Prüfe ob time verfügbar ist
    if ! command -v time &> /dev/null; then
        echo -e "${RED}Fehler: 'time' Kommando nicht gefunden!${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}✓ Alle Abhängigkeiten gefunden${NC}"
    echo ""
}

# Parameter prüfen
if [ $# -lt 1 ]; then
    print_usage
    exit 1
fi

PROGRAM_PATH="$1"
RUNS="${2:-3}"

# Prüfe ob Programm existiert und ausführbar ist
if [ ! -f "$PROGRAM_PATH" ]; then
    echo -e "${RED}Fehler: Programm '$PROGRAM_PATH' nicht gefunden!${NC}"
    exit 1
fi

if [ ! -x "$PROGRAM_PATH" ]; then
    echo -e "${RED}Fehler: Programm '$PROGRAM_PATH' ist nicht ausführbar!${NC}"
    echo "Versuche: chmod +x $PROGRAM_PATH"
    exit 1
fi

print_header
check_dependencies

PROGRAM_NAME=$(basename "$PROGRAM_PATH")
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="benchmark_results_${TIMESTAMP}"
mkdir -p "$LOG_DIR"

echo -e "${BLUE}Programm:${NC} $PROGRAM_NAME"
echo -e "${BLUE}Pfad:${NC} $PROGRAM_PATH"
echo -e "${BLUE}Anzahl Läufe:${NC} $RUNS"
echo -e "${BLUE}Log Verzeichnis:${NC} $LOG_DIR"
echo ""

# Arrays für Ergebnisse
declare -a runtimes
declare -a memory_peaks
declare -a cpu_usage

echo -e "${YELLOW}Starte Benchmarks...${NC}"
echo ""

for ((i=1; i<=RUNS; i++)); do
    echo -e "${BLUE}Lauf $i/$RUNS:${NC}"
    
    # Temporäre Dateien für diesen Lauf
    TIME_LOG="$LOG_DIR/time_$i.log"
    MEMORY_LOG="$LOG_DIR/memory_$i.log"
    
    # Runtime und Memory mit /usr/bin/time (detaillierte Statistiken)
    echo "  → Messe Runtime und Memory..."
    echo "  → Messe Runtime und Memory..."
    # Lese Parameter aus input.txt
    read -r fasta_path int_param < input.txt

    # Führe Programm mit Argumenten aus und messe Zeit
    /usr/bin/time -l "$PROGRAM_PATH" "$fasta_path" "$int_param"> "$MEMORY_LOG" > "$LOG_DIR/output_$i.log" 2> "$TIME_LOG"
    
    # Extrahiere Werte aus time output
    runtime=$(grep "real" "$TIME_LOG" | awk '{print $1}' | sed 's/real//')
    max_memory=$(grep "maximum resident set size" "$TIME_LOG" | awk '{print $1}')
    cpu_usage=$(grep "percent CPU" "$TIME_LOG" | awk '{print $1}' | sed 's/%//')
    
    # Konvertiere Memory von Bytes zu MB
    max_memory_mb=$(echo "scale=2; $max_memory / 1024 / 1024" | bc -l 2>/dev/null || echo "N/A")
    
    echo "    Runtime: $runtime"
    echo "    Peak Memory: ${max_memory_mb} MB"
    echo "    CPU Usage: ${cpu_usage}%"
    echo ""
    
    # Speichere Werte für Statistiken
    runtimes[$i]="$runtime"
    memory_peaks[$i]="$max_memory_mb"
    cpu_usage[$i]="$cpu_usage"
    
    # Kurze Pause zwischen den Läufen
    sleep 1
done

# Berechne Statistiken
echo -e "${GREEN}Berechne Statistiken...${NC}"

# Runtime Statistiken
runtime_sum=0
runtime_count=0
for rt in "${runtimes[@]}"; do
    if [[ "$rt" =~ ^[0-9]+\.[0-9]+$ ]]; then
        runtime_sum=$(echo "$runtime_sum + $rt" | bc -l)
        ((runtime_count++))
    fi
done

if [ $runtime_count -gt 0 ]; then
    avg_runtime=$(echo "scale=3; $runtime_sum / $runtime_count" | bc -l)
else
    avg_runtime="N/A"
fi

# Memory Statistiken
memory_sum=0
memory_count=0
for mem in "${memory_peaks[@]}"; do
    if [[ "$mem" =~ ^[0-9]+\.?[0-9]*$ ]] && [ "$mem" != "N/A" ]; then
        memory_sum=$(echo "$memory_sum + $mem" | bc -l)
        ((memory_count++))
    fi
done

if [ $memory_count -gt 0 ]; then
    avg_memory=$(echo "scale=2; $memory_sum / $memory_count" | bc -l)
else
    avg_memory="N/A"
fi

# Erstelle Zusammenfassungsreport
SUMMARY_FILE="$LOG_DIR/benchmark_summary.txt"
cat > "$SUMMARY_FILE" << EOF
C++ Programm Benchmark Report
============================
Programm: $PROGRAM_NAME
Pfad: $PROGRAM_PATH
Datum: $(date)
System: $(uname -a)
Anzahl Läufe: $RUNS

Ergebnisse:
-----------
Durchschnittliche Runtime: $avg_runtime Sekunden
Durchschnittlicher Peak Memory: $avg_memory MB

Einzelne Läufe:
EOF

for ((i=1; i<=RUNS; i++)); do
    echo "Lauf $i: Runtime=${runtimes[$i]}s, Memory=${memory_peaks[$i]}MB, CPU=${cpu_usage[$i]}%" >> "$SUMMARY_FILE"
done

# Zeige finale Ergebnisse
echo ""
echo -e "${GREEN}================================${NC}"
echo -e "${GREEN}        BENCHMARK ERGEBNISSE     ${NC}"
echo -e "${GREEN}================================${NC}"
echo ""
echo -e "${BLUE}Programm:${NC} $PROGRAM_NAME"
echo -e "${BLUE}Durchschnittliche Runtime:${NC} $avg_runtime Sekunden"
echo -e "${BLUE}Durchschnittlicher Peak Memory:${NC} $avg_memory MB"
echo ""
echo -e "${YELLOW}Detaillierte Logs gespeichert in:${NC} $LOG_DIR"
echo -e "${YELLOW}Zusammenfassung:${NC} $SUMMARY_FILE"
echo ""

# Öffne Zusammenfassung optional
read -p "Möchtest du die detaillierte Zusammenfassung anzeigen? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    cat "$SUMMARY_FILE"
fi

echo -e "${GREEN}Benchmark abgeschlossen!${NC}"