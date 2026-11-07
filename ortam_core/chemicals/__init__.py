"""Shared chemical database (toxicity values and physical-chemical properties).

Part of the shared core so that every module (risk, air, vadose, vapor, ...)
can read the same versioned data without calling another module.
Fortran analogy: a data MODULE that every program USEs.
"""
