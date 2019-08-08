# Read Mapping

An implementation of McCreight's Algorithm for read mapping.

## Instructions

Compile mcc.c and run executable with desired settings for FASTA file, read file, and alphabet. As many reference sequences and read counts are quite large (~5,300,000 characters, 500,000 reads), please allow time.

## Example

```bash
foo@bar~$ gcc -mcc.c
foo@bar~$ ./a.out fasta/Peach_reference.fasta reads/Peach_simulated_reads.fasta alphabet/DNA_alphabet.txt 
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


