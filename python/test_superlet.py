import numpy as np
from superlet import superlet, gen_superlet_testdata
# to get sth for the eyes ;)
import matplotlib.pyplot as ppl


def test_superlet():

    fs = 1000  # sampling frequency
    A = 10  # amplitude
    # 20Hz, 40Hz and 60Hz time localized harmonics with freq. neighbors
    signal = A * gen_superlet_testdata(freqs=[20, 40, 60], fs=fs, eps=0)

    # frequencies of interest in Hz
    foi = np.linspace(1, 100, 100)

    spec = superlet(
        signal,
        samplerate=fs,
        freqs=foi,
        order_max=30,
        order_min=1,
        c_1=5,
        adaptive=True,
    )

    # amplitude scalogram
    ampls = np.abs(spec)

    # tf profiles exactly at the harmonic frequencies
    assert foi[19] == 20
    assert foi[39] == 40
    assert foi[59] == 60

    # superlets taper off after 150ms
    t_off = 150

    # 20Hz is 'on' at around 0.5 seconds
    # and 'off' after around 1.5 seconds
    assert ampls[19][500] > 0.8 * A
    assert np.all(ampls[19][1500 + t_off:] < 0.2 * A)

    # 40Hz is 'on' at around 1.7 seconds
    # and 'off' after around 2.2 seconds
    assert ampls[39][1700] > 0.8 * A
    assert np.all(ampls[39][2200 + t_off:] < 0.2 * A)

    # 60Hz is 'on' at around 1.7 seconds
    # and 'off' before 2.25 seconds
    assert ampls[59][2520] > 0.8 * A
    assert np.all(ampls[59][:2250 - t_off] < 0.2 * A)

    fig1, ax1 = ppl.subplots(figsize=(5, 3.5))
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Scalogram profiles')
    tvec = np.arange(len(signal)) / fs
    ax1.plot(tvec, ampls[19], label='20Hz')
    ax1.plot(tvec, ampls[39], label='40Hz')
    ax1.plot(tvec, ampls[59], label='60Hz')
    ax1.plot([0, tvec[-1]], [A, A], 'k--', label='true amplitude')
    ax1.set_ylim((-.1, 15))
    ax1.legend(ncol=2)
    fig1.tight_layout()

    fig2, (ax2, ax3) = ppl.subplots(2, 1,
                                    sharex=True,
                                    gridspec_kw={"height_ratios": [1, 3]},
                                    figsize=(6, 6))

    ax2.plot(np.arange(signal.size) / fs, signal, c='cornflowerblue')
    ax2.set_ylabel('signal (a.u.)')

    extent = [0, len(signal) / fs, foi[0], foi[-1]]
    im = ax3.imshow(ampls, cmap="magma", aspect="auto", extent=extent, origin='lower')

    ppl.colorbar(im, ax=ax3, orientation='horizontal',
                 shrink=0.7, pad=0.2, label='amplitude (a.u.)')

    ax3.plot([0, len(signal) / fs], [20, 20], "--", c='0.5')
    ax3.plot([0, len(signal) / fs], [40, 40], "--", c='0.5')
    ax3.plot([0, len(signal) / fs], [60, 60], "--", c='0.5')

    ax3.set_xlabel("time (s)")
    ax3.set_ylabel("frequency (Hz)")

    fig2.tight_layout()

    return ampls
