Wizard
======

The OpenPathSampling Wizard is a friendly interactive tool for setting up
path sampling simulations. In general, the Wizard is best for users who are
new to path sampling and to Python -- it is very easy to use, but it is also
quite limited (relative to the core OPS library) in the functionality it can
handle.


Interacting with the Wizard
---------------------------

The Wizard is a text-based interactive guide to setting up a simulation. It
will ask you a question, and you give it an answer based on the simulation
you want to set up.

Sometimes the questions are multiple choice, in which case you can either
respond with the number of the choice or the text given in the choice:

.. code::

    🧙 What will you use for the underlying engine?

          1. OpenMM
          2. Load existing from OPS file

    🧙 Please select an option: 1

Sometimes the questions simply require you to fill in some text, such as a
file name:

.. code::

    🧙 Where is the XML file for your OpenMM integrator? integrator.xml

In some cases, you can use Python expressions, which the Wizard will
immediately calculate for you:

.. code::

    🧙 What is the minimum allowed value for 'psi' in this volume? 100 * np.pi / 180

For more on the abilities and limitation of these expressions, see
:ref:`Expression Evaluation <expression_eval>`.

Getting help in the Wizard
--------------------------

If you don't know what to do at one of the Wizard's prompts, input a ``?``
and the Wizard may give you a little more information about the question it
is asking:

.. code::

    🧙 What will you use for the underlying engine?

          1. OpenMM
          2. Load existing from OPS file

    🧙 Please select an option: ?

    🧙 An engine describes how you'll do the actual dynamics.

Sometimes you can ask the Wizard for more specific help by adding some text
after the question mark. For example, if the Wizard gives you a set of
choices, and you'd like to know more about one of the, follow the question
mark with the number or name of the choice:

.. code::

    🧙 What will you use for the underlying engine?

          1. OpenMM
          2. Load existing from OPS file

    🧙 Please select an option: ?1

    🧙 OpenMM is an GPU-accelerated library for molecular dynamics. To use
       OpenMM in the OPS wizard, you'll need to provide a file with topology
       information (usually a PDB), as well as XML versions of your OpenMM
       integrator and system objects.


Wizard commands
---------------

The Wizard has several commands, which are indicated by input beginning with
``!``, to allow you to control aspects of the wizard's workflow.  In
particular, the following commands can be invoked:

* ``!quit``: Quit the Wizard, with an option to save objects created so far.
* ``!!quit``: Immediately exit the Wizard (no saving of results).
* ``!restart``: Restart the current step in the Wizard. Useful if you
  mistype in a prompt, so that you can only restart this step of the wizard.


Name limitations in the Wizard
------------------------------

When using the Wizard, we recommend that names of objects not begin with
``?``, ``!``, and not be a numeral. This is because the Wizard frequently
presents users with a list of options which can be selected either by number
or by typing the exact text provided. Names that are numbers can cause
confusion between the numeric entry and the name, and names that begin with
``?`` or ``!`` would invoke the help/command systems instead of loading the
object. 
