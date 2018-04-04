# Accelerated HaplotypeCaller DNA Analysis Application on FPGAs Using CAPI



DNA carries all the information needed to define the genetic structure of an individual.
The large amounts of information stored in the DNA makes it costly to store and to process.
In this work, a fast and efficient implementation of a Field Programmable Gate
Array (FPGA) based, streaming multicore architecture for accelerating variant calling
algorithms will be designed. We focused on the HaplotypeCaller which is the variant
calling software part of the Genome Analysis Toolkit (GATK), whiich is one of the
most widely used DNA analysis tools in the field. The most time consuming part of
the HaplotypeCaller is the PairHMM algorithm. PairHMM is a probabilistic algorithm
that executes pairwise alignment of two sequences. Starting from an existing single
core PairHMM accelerator design, which was implemented using the POWER8 Coherent
Accelerator Processor Interface (CAPI), we designed three extra cores of the
PairHMM algorithm that can work independently and increase the performance of the
overall system. The new accelerator achieves a 2.2x speedup in comparison with the
single core. The accelerator is integrated with the HaplotypeCaller and uses CAPI to
access a shared processor-accelerator memory for direct communication. A JNI call is
used to send the memory addresses to the accelerator, which reduces the communication
overhead between the HaplotypeCaller and the accelerator. Results show that
the application is not able to saturate the accelerator with data, resulting in accelerator
idle time and under-utilization. The accelerator presented idle time for the 26% of the
different data sets that were used. It would be beneficial to implement a faster version
of the HaplotypeCaller which can load data sets to the accelerator as current data sets
are being processed.
