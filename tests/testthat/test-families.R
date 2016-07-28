context("families")

test_that("calculation of family probability matrices works", {
  kcs20.setoffamilies=by(names(kcs20), kcs20[,'type'], as.character)


  expect_equal_to_reference(kcs20.prob <- create_probab_sv_knowing_fam(kcs20.setoffamilies),
                            "testdata/kcs20.prob.rds")
  expect_equal(create_probab_sv_knowing_fam(kcs20.setoffamilies, db=kcs20), kcs20.prob)
})
