import yaml

meta = """\
package:
  name: {name}{dev}
  version: "{version}"

source:
#  git_url: ../../.git
  path: ../../

build:
  preserve_egg_dir: True
  number: 0

requirements:
  build:
    - python
    - setuptools
    - pyyaml

  run:
{requires}

test:
  requires:
    - nose
    - nose-timer
    - python-coveralls
    - msmbuilder
    - pyemma
    - ipynbtest

  imports:
    - openpathsampling

about:
  home: {home}
  license: {license}
  summary: {description}
"""


def main():
    # load settings from setup.py, easier to maintain
    with open('../../setup.yaml') as f:
        yaml_string = '\n'.join(f.readlines())

    prefs = yaml.load(yaml_string)

    print 'writing meta.yaml file in devtools/conda-recipe'

    tag = '    - '

    with open('../conda-recipe/meta.yaml', mode='w') as f:
        f.write(meta.format(
            name=prefs['name'],
            version=prefs['version'],
            dev='-dev',
            requires=tag + ('\n' + tag).join(prefs['requires']),
            home=prefs['download_url'],
            license=prefs['license'],
            description="'" + prefs['description'] + "'"
        ))


if __name__ == '__main__':
    main()
