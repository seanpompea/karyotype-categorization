create table dm_mpn.dbo.karyotype_classifications (
  mrn nvarchar(max)  null,
  result_date nvarchar(max)  null,
  karyotype nvarchar(max)  null,
  aml_swog nvarchar(max)  null,
  aml_calgb nvarchar(max)  null,
  aml_mrc_1998 nvarchar(max)  null,
  aml_mrc_2010 nvarchar(max)  null,
  mds_ipss nvarchar(max)  null,
  mds_5_group nvarchar(max)  null,
  mf_dipss_plus nvarchar(max)  null,
  is_normal nvarchar(max) null,
  abn_count nvarchar(max) null,
  full_numer_abn_tally nvarchar(max) null,
  trailing_mar_and_r_abn_tally nvarchar(max) null
);

