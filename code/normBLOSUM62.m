function B0 = normBLOSUM62()
    B0=blosum62;
    B0=B0(1:20,1:20);
    B0=B0+(abs(min(min((B0))))+1);
    B0=B0./max(max(B0));