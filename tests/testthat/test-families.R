context("families")

test_that("calculation of family probability matrices works", {
  kcs20.setoffamilies=by(names(kcs20), kcs20[,'type'], as.character)

  # Melina had the input neurons as the names of the vectors, which is not so
  # generally useful
  correct_families2=sapply(correct_families, names, simplify = F)
  expect_equal(create_probab_sv_knowing_fam(correct_families2[2:4]),
               probability_sv_knowing_family[,2:4], tol=1e-3)

  kcs20.prob <- create_probab_sv_knowing_fam(kcs20.setoffamilies)
  expect_equal_to_reference(kcs20.prob, "testdata/kcs20.prob.rds")
  expect_equal(create_probab_sv_knowing_fam(kcs20.setoffamilies, db=kcs20), kcs20.prob)
})
