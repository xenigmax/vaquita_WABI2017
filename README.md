This repositorty contains all the scripts that I used for benchmakring [Vaquita](https://github.com/seqan/vaquita).
The result was presented in [WABI 2017](http://www.acm-bcb.org/WABI/2017/index.php) in Boston. 

Prerequeistics
-----------------
* Variant callers

| Name         | Reference           |  Version  |
| -------------| ------------------- | --------- |
| Vaquita      | [source repository](https://github.com/seqan/vaquita) | WABI2017  |
| LumpyExpress | Layer et al., 2014  | 0.2.13    |
| Delly2       | Rausch et al., 2012 | 0.7.3     |
| CREST        | Wang et al., 2011   | 1.0       |
| Pindel       | Ye et al., 2009     | 0.2.5b8   |
| GASVPro      | Sindi et al., 2012  | Oct1_2013 |

You can use [wrapper scripts](https://github.com/xenigmax/vaquita_WABI2017/tree/master/wrapper) for several programs for your conveniecne.

* Simulators

| Name  | Reference           |  Version                   |
| ------| ------------------- | -------------------------- |
| Mason | Holtgrewe, 2010     | 2.0.5                      |
| ART   | Huang et al., 2012  | MountRainer-2016.06.05     |

You can refer my commands for generating simulation datasets in `simulations.txt`.

* Real datasets

| Name  | Link or Accession #         |
| ------| ------------------- |
| Genome in a bottle | ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed |
| Illumina Platinum Genomes   | ERR194147, ERR194160, and ERR194161 |


Configuration
-----------------
Please open `config.py` in this repository, and change settings accordingly.
For example, you may have to change binary files for variant callers in [here](https://github.com/xenigmax/vaquita_WABI2017/blob/master/config.py#L19-L26).

Run & Reporting
-----------------
You can refer my command in `commands.txt`.

Contact
-----------------
Jongkyu Kim (vaquita_WABI@jongkyu.kim)
