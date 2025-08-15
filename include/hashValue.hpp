#pragma once


struct hashValue
{   //=================================================================
    //KMER steht im komprimierten Suffix Array
    int cSAindex; // Index im komprimierten Suffix Array
    unsigned long occurences; // Anzahl der Vorkommen
    unsigned long lcp_interval_index; // Index des LCP Intervalls
    //=================================================================
    //REFERENZ auf anderes Pattern, von dem dieses abhängt
    int shift; // Verschiebung im Text (wenn -1 dann gibt keine Verschriebung d.h. nicht von anderem Pattern abhängig
    unsigned long refOccurrences;
    unsigned traceback_key; // Referenz auf das andere Pattern (kann nullptr sein, wenn es kein anderes Pattern gibt)
    bool processed;
               // sondern direkt in SA)

    
    hashValue()
    {
        cSAindex = -1; // noch nicht initialisiert
        occurences = 0; // noch nicht initialisiert
        shift = -1;  // noch nicht initialisiert
        refOccurrences = 0; // noch nicht 
        traceback_key = 0; // noch nicht initialisiert
        processed = false; // noch nicht initialisiert
    }

    void setValue (unsigned long cSAindex, unsigned long occurences)
    {
        this->cSAindex = cSAindex; 
        this->occurences += occurences;
        this -> processed = true;
    }

    void setReferenceValue(int shift, unsigned long refOccurrences, unsigned traceback_key, bool processed)
    {
        this->shift = shift;
        this->refOccurrences = refOccurrences;
        this->traceback_key = traceback_key;
        this->processed = processed;
    }

    
};