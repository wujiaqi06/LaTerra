\name{lt_rate_common_domain}
\alias{lt_rate_common_domain}
\alias{validate_lt_rate_common_domain}
\alias{lt_compare_rate_pair}
\alias{lt_c3_c2_compare}
\title{Exact keyed common domains and C3-vs-C2 sensitivity comparison}
\usage{
lt_rate_common_domain(...)

validate_lt_rate_common_domain(x)

lt_compare_rate_pair(first, second, positive_ratio = FALSE,
  comparison = paste(first$representation, second$representation, sep = "_vs_"))

lt_c3_c2_compare(c3_fit, c2_fit)
}
\arguments{
  \item{...}{At least two \code{lt_rate} objects, or one list of them.}
  \item{x}{An \code{lt_rate_common_domain}.}
  \item{first, second}{Rate objects sharing exact keyed source/coordinate/Layer-2 identity.}
  \item{positive_ratio}{Whether positive-ratio summaries are meaningful.}
  \item{comparison}{Stable comparison label.}
  \item{c3_fit}{Certified production C3 fit.}
  \item{c2_fit}{Certified required-sensitivity C2 fit.}
}
\description{
Constructs exact keyed availability intersections and descriptive paired
summaries. Direct C3-vs-C2 comparison uses the common positive domain only.
}
\details{
Full-domain results with different availability are labelled
\code{UNEQUAL_DOMAIN_DESCRIPTIVE}. Domain loss is not estimator-performance
evidence. C3-vs-C2 divergence is the sensitivity result; no winner or biological
validation is declared. The first-class common-domain object stores compact
scientific axes and \code{common_linear}; validation reconstructs its common
mask and keyed identity. Mutable display labels are outside that identity.
}
\value{
An \code{lt_rate_common_domain}, a paired-summary data frame, or an
\code{lt_c3_c2_comparison}; validation returns its input invisibly.
}
