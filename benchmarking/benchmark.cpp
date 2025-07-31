#include <benchmark/benchmark.h>
#include <sys/resource.h>
#include "compressV2.hpp"
#include "fastaParser.hpp"
#include "suffix_array.hpp"
#include <vector>
#include <string>

// Globale Variable für die Dateipfade
static std::vector<std::string> g_test_files;

// Hilfsfunktion zur Speichermessung
long getPeakMemoryUsageKB() {
    struct rusage usage;
    // RUSAGE_SELF fragt den Speicher für den aktuellen Prozess ab
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // Auf macOS ist ru_maxrss in Bytes, auf Linux in Kilobytes.
        // Wir müssen das zur Portabilität berücksichtigen.
#if defined(__APPLE__) && defined(__MACH__)
        return usage.ru_maxrss / 1024; // Umrechnung von Bytes nach KB
#else
        return usage.ru_maxrss; // Ist bereits in KB
#endif
    }
    return 0;
}

// Die korrigierte Benchmark-Funktion
void BM_uncompressedSA(benchmark::State& state) {
    // 1. Einmaliges Setup (wird nicht gemessen)
    long file_index = state.range(0);
    std::string fastaData = parseFasta(g_test_files[file_index]);

    // Misst den Speicher vor Beginn der Messschleife
    long mem_before = getPeakMemoryUsageKB();

    // 2. Die eigentliche Messschleife (misst die Zeit)
    for (auto _ : state) {
        // Das ist die Operation, deren Zeit wir messen wollen.
        SuffixArray SA(fastaData);

        // Verhindert, dass der Compiler die Erstellung des Objekts wegoptimiert.
        // benchmark::DoNotOptimize(compressedSA);
    }

    // 3. Setzen der Counter (nach der Schleife!)
    long mem_after = getPeakMemoryUsageKB();
    // Wir berichten den absoluten Peak-Speicher nach der Ausführung.
    state.counters["Peak RSS (KB)"] = mem_after;
    // Optional: Berichte auch die Differenz, falls das für dich interessant ist.
    state.counters["Delta RSS (KB)"] = mem_after - mem_before;

    SuffixArray SA(fastaData);
    state.counters["Exact Memory (Byte)"] = SA.memoryUsageBytes(); // Umrechnung in KB
}



// Die korrigierte main()-Funktion
int main(int argc, char** argv) {
    // Wir müssen die Argumente von Google Benchmark trennen.
    // Alle Argumente, die nicht mit '-' beginnen, sind unsere Dateipfade.
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            g_test_files.push_back(argv[i]);
            // "Entferne" das Argument, damit Google Benchmark es nicht sieht.
            // Dies ist ein gängiger Trick: überschreibe es mit dem letzten und kürze die Liste.
            argv[i] = argv[argc - 1];
            argc--;
            i--; // Überprüfe das neue Argument an dieser Position erneut
        }
    }

    if (g_test_files.empty()) {
        fprintf(stderr, "Fehler: Keine Eingabedateien angegeben.\n");
        return 1;
    }

    // Dynamische Registrierung für jede gefundene Datei
    for (int i = 0; i < g_test_files.size(); ++i) {
        benchmark::RegisterBenchmark("BM_uncompressedSA", &BM_uncompressedSA)->Arg(i);
    }

    // Google Benchmark initialisieren und ausführen
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    
    return 0;
}