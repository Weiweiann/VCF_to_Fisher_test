# VCF_to_Fisher_test

Change the "GT" value into genotype in .vcf files.

The file name of Controls and Patients can be given in ```.yaml``` file. As the
example shows below.

```
Controls:
- absolute_path_of_control1.vcf
- absolute_path_of_control2.vcf
- absolute_path_of_control3.vcf

Patients:
- absolute_path_of_patients1.vcf
- absolute_path_of_patients2.vcf
- absolute_path_of_patients3.vcf
```

The fisher exact test will follow the group definition and caculating the
p-value and odds ratio, and save the result as ``Fisher_test_output.csv``




### Notice
This version doesnot handle the vcf with indel but only snps.

