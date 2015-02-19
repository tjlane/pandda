from scipy.optimize import fsolve
import numpy

def variable_sigma_d_log_likelihood(est_sigma, est_mu, obs_vals, obs_error):
    """Calculate the value of the differentiated log likelihood for the values of mu, sigma"""
    term1 = (obs_vals - est_mu)**2 / ((est_sigma**2 + obs_error**2)**2)
    term2 = 1 / (est_sigma**2 + obs_error**2)
    return numpy.sum(term1) - numpy.sum(term2)

def estimate_true_underlying_mean_and_sd(obs_vals, obs_error):
    """Given a  set of observations `obs_vals` with estimated errors `obs_sigma`, estimate the mean and sd of the underlying distribution"""

    obs_vals  = numpy.array(obs_vals)
    obs_error = numpy.array(obs_error)

    # First calculate the mean of the sample - this is a good estimator of the true mean
    est_mu = numpy.mean(obs_vals)

    # Initial guess
    est_sigma = 1

    return fsolve(func=variable_sigma_d_log_likelihood, x0=est_sigma, args=(est_mu, obs_vals, obs_error))


if __name__ == '__main__':

    # True values we are trying to estimate
    true_mean = 10
    true_sd   = 5

    # Number of observations
    num_obs = 500

    guesses  = []

    for attempt in xrange(200):

        print('============================>')

        # Sample the original distribution
        true_vals = true_mean + true_sd*numpy.random.randn(num_obs)
        print('MEAN OF TRUE VALS: {!s} ({!s})'.format(round(numpy.mean(true_vals),3), true_mean))
        print(' STD OF TRUE VALS: {!s} ({!s})'.format(round(numpy.std(true_vals),3), true_sd))

        # Create a random selection of sigmas for the different observations
        obs_error = numpy.abs(2 + 1*numpy.random.randn(num_obs))
        print('MEAN OF OBS ERROR: {!s} ({!s})'.format(round(numpy.mean(obs_error),3), '2ish'))
        print(' STD OF OBS ERROR: {!s} ({!s})'.format(round(numpy.std(obs_error),3), '1ish'))

        # Noise to be added to the true observations
        obs_noise = numpy.random.randn(num_obs)
        print('MEAN OF OBS NOISE: {!s} ({!s})'.format(round(numpy.mean(obs_noise),3), 0))
        print(' STD OF OBS NOISE: {!s} ({!s})'.format(round(numpy.std(obs_noise),3), 1))

        # Create fake data!
        obs_vals = true_vals + obs_error * obs_noise
        print(' MEAN OF OBS VALS: {!s}'.format(round(numpy.mean(obs_vals),3)))
        print('  STD OF OBS VALS: {!s}'.format(round(numpy.std(obs_vals),3)))

        out = estimate_true_underlying_mean_and_sd(obs_vals=obs_vals, obs_error=obs_error)

        print(' ESTIMATED VALUES: {!s}'.format(round(out,3)))

        guesses.append(out)

