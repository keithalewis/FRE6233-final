# FRE6233 Final Exam

You have forked https://github.com/keithalewis/FRE6233-final 
to your GitHub repository and cloned that to your local machine.
Submit your solutions by pushing to your fork on GitHub and
issuing a pull request to have your work graded.

Written Exam

Write answers in this file or on paper.

Let (B_t) be standard Brownian motion.
The Bachelier model of stock price is F_t = f(1 + sigma B_t),
where f is the forward and sigma is the volatility

1. What is E[F_t] and Var(F_t)?

Solution:

Consider a European option that pays (F_T)^2 at time T.

2. Derive a formula for V_t = E_t[F_T^2] in terms of F_t.

Hint. Use B_t^2 - t is a martingale.

Solution:

Let Z be a standard normal random variable and N(z) = P(Z <= z) be its cumulative distribution function.

3. Show E[Z 1(Z <= a)] = -n(a) where n(z) = N'(z) is the standard normal density function.

Solution:

4. Prove the forward value of a European put option with strike k and expiration T is

    p = E[max{k - F_T, 0}] = (k - f) N(x) + f sigma sqrt(t) n(x)

where x = (k - f)/(f sigma sqrt(t)).

Hint. B_t/sqrt(t) is standard normal.

Solution:

5. Derive a closed form solution for the Bachelier delta dp/df.

Hint. Use dg(F_T)/df = g'(F_T) dF_T/df where g(x) = max{k - x, 0}.

Solution:

Coding Exam

Follow the instructions in xll_template.cpp.
