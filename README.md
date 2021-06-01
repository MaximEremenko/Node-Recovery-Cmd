# Node Recovery Cmd

WIDE COLUMN ERASURE CODE FOR LARGE CAPACITY DISTRIBUTED STORAGE SYSTEM.

The main goal of this project is to propose codes that optimize node repair in terms of reconstruction rate that is amount of downloaded data during reconstruction divided by amount of surviving data. Proposed codes must ensure high performance in terms of encoding/decoding/reconstruction execution time. To ensure scalability proposed codes must support dynamic stripe scaling, such as capacity expansion and recombination. During capacity expansion and recombination, the algorithms keep the parity nodes number unchanged. When adding or deleting a data node, the algorithms must guarantee that the network overhead is optimal.
